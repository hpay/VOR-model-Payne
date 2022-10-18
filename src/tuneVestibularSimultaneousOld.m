function [K0new, K1new, K2new] = tuneVestibularSimultaneousOld(...
    K0_0, K1_0, K2_0, B, I, S, ...
    tts, head, target, E, P, conds, light, sines, R, RL_data)

% PARAMETERS
debug_on = 0;
% scale_regularization = 1; % scale the total regularization penalty
dt = K1_0.dt;
mask_learn = B.learn_mask_vest_eye_pc;
mask_tune = B.tune_mask_vest_eye_pc;
n_total = length(K1_0.Kb_eye_pc);
n_final = nnz(mask_tune);
% END PARAMETERS

% Goal gain and phase of eye
E_gain_goal1 = RL_data.gains1/RL_data.scaleJR2SL;
E_phase_goal1 = RL_data.phases1;
E_gain_goal2 = RL_data.gains2/RL_data.scaleJR2SL;
E_phase_goal2 = RL_data.phases2;
E_gain_goal0 = RL_data.gains0/RL_data.scaleJR2SL;
E_phase_goal0 = RL_data.phases0;

% Get the goal PC baseline gain and phase for 0.5 Hz (close to steady state)
[pc_amp, pc_phase] = fitsine(dt, P{1}, .5);

vord_step_ind = strcmp(conds,'500ms_dark');
pc_gain_ss = I.pc_gain_ss;
% pc_gain_ss = (nanmean(P{vord_step_ind}(400:600))/nanmean(head{vord_step_ind}(400:600)));
trans_mask = tts{vord_step_ind}>0 & tts{vord_step_ind}< I.T_step_trans;
temp1 = P{vord_step_ind}(trans_mask);
pc_trans_init = temp1(find(abs(temp1)==max(abs(temp1)),1))/max(head{vord_step_ind});

eye_ss = [ E_gain_goal0(1) E_gain_goal1(1) E_gain_goal2(1)];

% Combining SL and Watanabe
pc_ss = pc_gain_ss + 1/RL_data.scaleJR2SL*...
    I.pc_ss.*(1-I.eye_ss); % change in sp/s per change in deg/s, * change in gain


% Create a PC goal gain and phase
P_gain_goal1 = NaN(size(RL_data.freqs));
P_phase_goal1 = NaN(size(RL_data.freqs));
P_gain_goal1(1) = abs(pc_gain_ss);
P_phase_goal1(1) = pc_phase;

P_gain_goal0 = NaN(size(RL_data.freqs));
P_phase_goal0 = NaN(size(RL_data.freqs));
P_gain_goal0(1) = abs(pc_ss(1));
P_phase_goal0(1) = 11; % deg, Watanabe 1985

P_gain_goal2 = NaN(size(RL_data.freqs));
P_phase_goal2 = NaN(size(RL_data.freqs));
P_gain_goal2(1) = abs(pc_ss(3));
P_phase_goal2(1) = 172; % deg, Watanabe 1985


% Using Lisberger data for the transients (don't weight too much, likely to
% ohigh because it measures the peak rate in 1st 50 ms)
pc_trans =  pc_trans_init + 1/RL_data.scaleJR2SL*[1.62  0 -1.01];% TODO: check the scaling

% Original free parameters (3 x n)
p_init = [K0_0.Kb_eye_pc(mask_tune)...
    K1_0.Kb_eye_pc(mask_tune)...
    K2_0.Kb_eye_pc(mask_tune)];

% Find the initial net weight of kEP to enforce that it doesn't change
sumEP = sum(K1_0.EP);
temp_tune = find(mask_tune);
temp_learn = find(mask_learn);
mask_fixEP = ~ismember(temp_tune, temp_learn);

% Set lower and upper bounds
lb = -inf(n_total, 3);
ub = inf(n_total, 3);



% Add a constraint so that baseline kEP is positive
if I.restrict_EP_pos
    lb(B.bEP_mask,:) = 0;
end

% Mask out the lower and upper bounds
lb = lb(mask_tune,:);
ub = ub(mask_tune,:);

% Plot progress
% figure; hs = subplot(1,1,1);
% set(hs, 'ColorOrder', parula(10), 'NextPlot', 'replacechildren');

legend_str = [];

% Initialize these so they are global variables for plotting
err_plot = [];
err_out_curr = [];

% Set up A and B constraints to constrain the direction of plasticity at kPH

% Constrain the NET direction of plasticity if specified: A_ineq*X<b_ineq
A_ineq = [ ...
    -I.fixPH*(B.bPH_mask)   I.fixPH*(B.bPH_mask)    zeros(1,n_total);...
    zeros(1,n_total)        I.fixPH*(B.bPH_mask)    -I.fixPH*(B.bPH_mask)];

% Only include cols for the weights that are being tuned
A_ineq = A_ineq(:,...
    [mask_tune(:); mask_tune(:); mask_tune(:)]);

% Mask out the rows that are not constrainted (all zero)
A_ineq = A_ineq(any(A_ineq,2),:);

B_ineq = zeros(size(A_ineq,1),1);   % B is all zeros


% METHOD 2: fmincon
options = optimoptions('fmincon','MaxIter',1000,'FunctionTolerance',10,...
    'StepTolerance',1e-7,'MaxFunctionEvaluations', 1e8,...
    'Display','iter','UseParallel',true,'PlotFcn','optimplotfval');

[p_out, fval, exitflag, output] = fmincon(@(x)myMinFun(x, I), p_init,A_ineq,B_ineq,...
    [],[], lb,ub, [], options);
disp(output.message)
p0 = p_out(:,1);
p1 = p_out(:,2);
p2 = p_out(:,3);

% Enfore kEP isn't changing with learning - force total weight of kEP to stay the same
K1new_init = updateKp(K1_0, I, p1, mask_tune);        % Temporary variable
sumEP_curr = sum(K1new_init.EP);
p1(mask_fixEP) =  p1(mask_fixEP) * sumEP/sumEP_curr;

% Enfore kEP isn't changing with learning
p0(mask_fixEP) = p1(mask_fixEP);
p2(mask_fixEP) = p1(mask_fixEP);

K0new = updateKp(K1_0, I, p0, mask_tune);
K1new = updateKp(K1_0, I, p1, mask_tune);
K2new = updateKp(K1_0, I, p2, mask_tune);

    function err_out = myMinFun(p, I)
        
        p0_curr = p(:,1);
        p1_curr = p(:,2);
        p2_curr = p(:,3);
        
        % Enfore kEP isn't changing with learning: Force total weight of kEP to stay the same
        K1_0_init = updateKp(K1_0, I, p1_curr, mask_tune);        % Temporary variable
        sumEP_curr = sum(K1_0_init.EP);
        p1_curr(mask_fixEP) =  p1_curr(mask_fixEP) * sumEP/sumEP_curr;
        
        % Enfore kEP isn't changing with learning
        p0_curr(mask_fixEP) = p1_curr(mask_fixEP);
        p2_curr(mask_fixEP) = p1_curr(mask_fixEP);
        
        % Update structure K with new variable weights
        K0_0 = updateKp(K0_0, I, p0_curr, mask_tune);
        K1_0 = updateKp(K1_0, I, p1_curr, mask_tune);
        K2_0 = updateKp(K2_0, I, p2_curr, mask_tune);
        
        % Get impulse if needed
        if I.impulse_or_closedloop
            KF0 = getImpulse(K0_0, I, 1);       % 1 = dark filters only, faster
            KF1 = getImpulse(K1_0, I, 1);       % 1 = dark filters only, faster
            KF2 = getImpulse(K2_0, I, 1);       % 1 = dark filters only, faster
        else
            KF0 = K0_0;
            KF1 = K1_0;
            KF2 = K2_0;
        end
        
        % Get errors
        [err_eye_all0, err_pc_all0, err_reg0] = getError(I, p0_curr, p1_curr, KF0,KF1,0);
        [err_eye_all1, err_pc_all1, err_reg1] = getError(I, p1_curr, p1_curr, KF1,KF1,1);
        [err_eye_all2, err_pc_all2, err_reg2] = getError(I, p2_curr, p1_curr, KF2,KF1,2);
        
        % Sum error for each learning condition
        err_out0 = [nansum(err_eye_all0)  nansum(err_pc_all0) I.scale_learn*err_reg0];
        err_out1 = [nansum(err_eye_all1)  nansum(err_pc_all1)  err_reg1];
        err_out2 = [nansum(err_eye_all2)  nansum(err_pc_all2)  I.scale_learn*err_reg2];
        
        % Error for deviation from original EP weight
        err_EP  = I.scale_const_EP * (sumEP_curr - sumEP).^2;
        
        % Total error
        err_out_curr = [err_out0  err_out1  err_out2  err_EP];
        legend_str = {'E0','P0','reg0','E1','P1','reg1','E2','P2','reg2','scaleEP'};
        err_out = sum(err_out_curr);
    end


    function [err_eye, err_pc, err_reg] = getError(I, pX_curr, p1_curr, KF_X, KF_1, type)
        
    
%         plotBaselineResults(KF_1, I, conds, tts, head, target, E, P, light, sines, [1 2 3 4 23 24 25]);
 
        % Get the current frequency response
        [Ehat_freq, Phat_freq, E_gain, E_phase, P_gain, P_phase] = getFreq(KF_X, I, RL_data.freqs);
        [Ehat_freq1, Phat_freq1, E_gain1, E_phase1, P_gain1, P_phase1] = getFreq(KF_1, I, RL_data.freqs);
        E_phase1 = wrapTo180(E_phase1-180);
        E_phase = wrapTo180(E_phase-180);
        
        % Frequency error: difference in gain and phase from desired
        if type==1
            err_eye_freq = ((log(E_gain)- log(E_gain_goal1)).^2 + ...
                deg2rad(wrapTo180(E_phase - E_phase_goal1)).^2);
            err_pc_freq = ((log(P_gain)- log(P_gain_goal1)).^2 + ...
                deg2rad(wrapTo180(P_phase - P_phase_goal1)).^2);
            
        elseif type==0
            err_eye_freq = ( ( (log(E_gain)-log(E_gain1)) - (log(E_gain_goal0)-log(E_gain_goal1) ) ).^2 + ...
                deg2rad(wrapTo180(  (E_phase-E_phase1) - (E_phase_goal0-E_phase_goal1) )).^2 );
            err_pc_freq = ((log(P_gain)- log(P_gain_goal0)).^2 + ...
                deg2rad(wrapTo180(P_phase - P_phase_goal0)).^2);
            
        elseif type==2
            err_eye_freq = ( ( (log(E_gain)-log(E_gain1)) - (log(E_gain_goal2)-log(E_gain_goal1) ) ).^2 + ...
                deg2rad(wrapTo180(  (E_phase-E_phase1) - (E_phase_goal2-E_phase_goal1) )).^2 );
            err_pc_freq = ((log(P_gain)- log(P_gain_goal2)).^2 + ...
                deg2rad(wrapTo180(P_phase - P_phase_goal2)).^2);
            
        end
        
        
        % Baseline error: JR data in the dark
        SSE = @(xhat, x) cellfun(@(xhat, x) nansum((xhat - x).^2), xhat, x);
        if type == 1
            mask_JR = ~light & (1:length(light))'<=I.nConds_JR;
            [Ehat_JR, Phat_JR] = getModelJR(KF_X, I, head(mask_JR), target(mask_JR), light(mask_JR), sines(mask_JR));
            
            err_eye_JR = SSE(Ehat_JR, E(mask_JR));
            err_pc_JR = SSE(Phat_JR, P(mask_JR));
        else
            % Specific to learning: only look at Ramachandran & Lisberger frequency data
            err_eye_JR = [];
            err_pc_JR = [];
        end
        
        % STEP
        step_ind = strcmp(conds,'step');
        tt_step = tts{step_ind};
        head_step = head{step_ind}/max(head{step_ind});
        mask_step_ss =  tt_step > I.T_step_ss;  % Start comparing 500 ms in
        mask_step_trans =  tt_step > 0 & tt_step < I.T_step_trans; % Transient is the peak prior to 50 ms - increase to
        
        [Ehat_step, Phat_step] = getModelJR(KF_X, I, {head_step}, {zeros(size(tt_step))}, 0, 0);
        err_eye_step = sum( ( -Ehat_step(mask_step_ss) - eye_ss(type+1)).^2 );
        err_pc_step  = sum( ( Phat_step(mask_step_ss) - pc_ss(type+1)).^2 );
        
        % Get the error in the transient
        if type ~= 1
            temp = Phat_step(mask_step_trans);
            pc_trans_curr = temp(find(abs(temp)==max(abs(temp)),1));
            err_pc_step_trans = (pc_trans_curr - pc_trans(type+1)).^2; % Transient = peak firing in 1st 50 ms after step onset
        else
            % Don't penalize the transient if we are before learning
            err_pc_step_trans  = 0;
        end
        weight_freqs = 100000;
        weight_step = 10000; % changes from 1000
        weight_transient = 10000; % Don't weight this too much - max measurement of Lisberger 1994 II. may exagerate transient in the data
        
        % TODO: check scale
        err_eye = [err_eye_JR;...
            weight_freqs.*err_eye_freq;... % D
            weight_step*err_eye_step];              % Step steady state
        err_pc = [err_pc_JR; ...
            weight_freqs.*err_pc_freq; ...
            weight_step*err_pc_step;...
            weight_transient*err_pc_step_trans]; % TODO check this!
        
        
        % REGULARIZATION ERROR
        % Create combine regularization matrices for both eye and PC params
        Tall0  = blkdiag(R.T0_eye, R.T0_pc);
        Tall1  = blkdiag(R.T1_eye, R.T1_pc);
        Tall2  = blkdiag(R.T2_eye, R.T2_pc);
        
        % Mask out only the params we care about now
        Tall0_mask = Tall0(mask_tune, mask_tune);
        Tall1_mask = Tall1(mask_tune, mask_tune);
        Tall2_mask = Tall2(mask_tune, mask_tune);
        
        if type == 1
            err_reg_all = I.reg_lambda0*(Tall0_mask*p1_curr).^2 +...
                I.reg_lambda1*(Tall1_mask*p1_curr).^2 + I.reg_lambda2*(Tall2_mask*p1_curr).^2;
        else
            % Penalize CHANGE in weights for learning
            err_reg_all = I.reg_lambda0*(Tall0_mask*(pX_curr-p1_curr)).^2 + ...
                I.reg_lambda1*(Tall1_mask*(pX_curr-p1_curr)).^2 + I.reg_lambda2*(Tall2_mask*(pX_curr-p1_curr)).^2;
        end
        err_reg = sum(err_reg_all);
        
    end

    function K_out = updateKp(K_in, I, p_new, mask_tune)
        Kb_eye_pc_new = K_in.Kb_eye_pc;
        Kb_eye_pc_new(mask_tune) = p_new;
        
        Kb_eye = Kb_eye_pc_new(1:length(K_in.Kb_eye));
        Kb_pc = Kb_eye_pc_new(length(K_in.Kb_eye)+1:end);
        
        K_out = updateK(B, I, S,  Kb_eye, Kb_pc);
    end



end

