function [K0new, K1new, K2new] = tuneVestibularSimultaneous(...
    K0_0, K1_0, K2_0, B, I, S, ...
    tts, head, target, E, P, conds, light, sines, R, RR_data)

% Mask free params
mask_learn = B.learn_mask_vest_eye_pc;
mask_tune = B.tune_mask_vest_eye_pc;
n_total = length(K1_0.Kb_eye_pc);

% Combine all JR's data in one structure
JR.conds = conds;
JR.tts = tts;
JR.head = head;
JR.target = target;
JR.E = E;
JR.P = P;
JR.light = light;
JR.sines = sines;

% Goal gain and phase of eye
RR.E_gain_goal1 = RR_data.gains1/RR_data.scaleJR2RR;
RR.E_phase_goal1 = RR_data.phases1;
RR.E_gain_goal2 = RR_data.gains2/RR_data.scaleJR2RR;
RR.E_phase_goal2 = RR_data.phases2;
RR.E_gain_goal0 = RR_data.gains0/RR_data.scaleJR2RR;
RR.E_phase_goal0 = RR_data.phases0;

% Get the goal PC baseline gain and phase for 0.5 Hz (close to steady state)
[pc_amp, pc_phase] = fitsine(K1_0.dt, P{1}, .5);

vord_step_ind = strcmp(conds,'500ms_dark');
pc_gain_ss = I.pc_gain_ss;
% pc_gain_ss = (nanmean(P{vord_step_ind}(400:600))/nanmean(head{vord_step_ind}(400:600)));
trans_mask = tts{vord_step_ind}>0 & tts{vord_step_ind}< I.T_step_trans;
temp1 = P{vord_step_ind}(trans_mask);
pc_trans_init = temp1(find(abs(temp1)==max(abs(temp1)),1))/max(head{vord_step_ind});

WL.eye_ss = [ RR.E_gain_goal0(1) RR.E_gain_goal1(1) RR.E_gain_goal2(1)];
%       eye_ss = [ E_gain_goal0(1) E_gain_goal1(1) E_gain_goal2(1)];

% Combining SL and Watanabe
WL.pc_ss = pc_gain_ss + 1/RR_data.scaleJR2RR*...
    I.pc_ss.*(1-I.eye_ss); % change in sp/s per change in deg/s, * change in gain


% Create a PC goal gain and phase
RR.P_gain_goal1 = NaN(size(RR_data.freqs));
RR.P_phase_goal1 = NaN(size(RR_data.freqs));
RR.P_gain_goal1(1) = abs(pc_gain_ss);
RR.P_phase_goal1(1) = pc_phase;

RR.P_gain_goal0 = NaN(size(RR_data.freqs));
RR.P_phase_goal0 = NaN(size(RR_data.freqs));
RR.P_gain_goal0(1) = abs(WL.pc_ss(1));
RR.P_phase_goal0(1) = 11; % deg, Watanabe 1985

RR.P_gain_goal2 = NaN(size(RR_data.freqs));
RR.P_phase_goal2 = NaN(size(RR_data.freqs));
RR.P_gain_goal2(1) = abs(WL.pc_ss(3));
RR.P_phase_goal2(1) = 172; % deg, Watanabe 1985

RR.freqs = RR_data.freqs;
 
% Using Lisberger data for the transients
WL.pc_trans =  pc_trans_init + 1/RR_data.scaleJR2RR*[1.62  0 -1.01];% ***TODO: check the scaling

% Original free parameters (3 x n)
p_init = [K0_0.Kb_eye_pc(mask_tune)...
    K1_0.Kb_eye_pc(mask_tune)...
    K2_0.Kb_eye_pc(mask_tune)];

% Find the initial net weight of kEP to enforce that it doesn't change
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

% Set options for fmincon
options = optimoptions('fmincon','MaxIter',1000,'FunctionTolerance',10,...
    'StepTolerance',1e-7,'MaxFunctionEvaluations', 1e8,...
    'Display','iter','UseParallel',true,'PlotFcn','optimplotfval');

[p_out, fval, exitflag, output] = fmincon(@(x)myMinFun(x, B, I, S,...
    K1_0, K2_0, K0_0, mask_tune, mask_fixEP, JR, RR, WL, R), p_init,A_ineq,B_ineq,...
    [],[], lb,ub, [], options);
disp(output.message)
p0 = p_out(:,1);
p1 = p_out(:,2);
p2 = p_out(:,3);

%% TODO: check unnecessary??
% Enfore kEP isn't changing with learning - force total weight of kEP to stay the same
K1new_init = updateKp(K1_0, B, I, S, p1, mask_tune);        % Temporary variable
sumEP_final = sum(K1new_init.EP);
sumEP_orig = sum(K1_0.EP);
p1(mask_fixEP) =  p1(mask_fixEP) * sumEP_orig/sumEP_final;

% Enfore kEP isn't changing with learning
p0(mask_fixEP) = p1(mask_fixEP);
p2(mask_fixEP) = p1(mask_fixEP);

K0new = updateKp(K1_0, B, I, S, p0, mask_tune);
K1new = updateKp(K1_0, B, I, S, p1, mask_tune);
K2new = updateKp(K1_0, B, I, S, p2, mask_tune);

%%
    function err_out = myMinFun(p, B, I, S, K1_0, K2_0, K0_0, mask_tune, mask_fixEP, JR, RR, WL, R)
        
        p0_curr = p(:,1);
        p1_curr = p(:,2);
        p2_curr = p(:,3);
        
        % Enfore kEP isn't changing with learning: Force total weight of kEP to stay the same
        sumEP = sum(K1_0.EP);
        K1_0_init = updateKp(K1_0, B, I, S, p1_curr, mask_tune);        % Temporary variable
        sumEP_curr = sum(K1_0_init.EP);
        p1_curr(mask_fixEP) =  p1_curr(mask_fixEP) * sumEP/sumEP_curr;
        
        % Enfore kEP isn't changing with learning
        p0_curr(mask_fixEP) = p1_curr(mask_fixEP);
        p2_curr(mask_fixEP) = p1_curr(mask_fixEP);
        
        % Update structure K with new variable weights
        K0_0 = updateKp(K0_0, B, I, S, p0_curr, mask_tune);
        K1_0 = updateKp(K1_0, B, I, S, p1_curr, mask_tune);
        K2_0 = updateKp(K2_0, B, I, S, p2_curr, mask_tune);
        
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
        [err_eye_all0, err_pc_all0, err_reg0] = getError(I, p0_curr, p1_curr, KF0,KF1,0, mask_tune, JR, RR, WL, R);
        [err_eye_all1, err_pc_all1, err_reg1] = getError(I, p1_curr, p1_curr, KF1,KF1,1, mask_tune, JR, RR, WL, R);
        [err_eye_all2, err_pc_all2, err_reg2] = getError(I, p2_curr, p1_curr, KF2,KF1,2, mask_tune, JR, RR, WL, R);
        
        % Sum error for each learning condition
        err_out0 = [nansum(err_eye_all0)  nansum(err_pc_all0) I.scale_learn*err_reg0];
        err_out1 = [nansum(err_eye_all1)  nansum(err_pc_all1)  err_reg1];
        err_out2 = [nansum(err_eye_all2)  nansum(err_pc_all2)  I.scale_learn*err_reg2];
        
        % Error for deviation from original EP weight
        err_EP  = I.scale_const_EP * (sumEP_curr - sumEP).^2;
        
        % Total error
        err_out_curr = [err_out0  err_out1  err_out2  err_EP];
%         legend_str = {'E0','P0','reg0','E1','P1','reg1','E2','P2','reg2','scaleEP'};
        err_out = sum(err_out_curr);
    end


    function [err_eye, err_pc, err_reg] = getError(I, pX_curr, p1_curr, KF_X, KF_1, type, mask_tune,JR, RR, WL, R)
        
        % Get the current frequency response
        [Ehat_freq, Phat_freq, E_gain, E_phase, P_gain, P_phase] = getFreq(KF_X, I, RR.freqs);
        [Ehat_freq1, Phat_freq1, E_gain1, E_phase1, P_gain1, P_phase1] = getFreq(KF_1, I, RR.freqs);
        E_phase1 = wrapTo180(E_phase1-180);
        E_phase = wrapTo180(E_phase-180);
        
        % Frequency error: difference in gain and phase from desired
        if type==1
            err_eye_freq = ((log(E_gain)- log(RR.E_gain_goal1)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(E_phase - RR.E_phase_goal1)).^2);
            err_pc_freq = ((log(P_gain)- log(RR.P_gain_goal1)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(P_phase - RR.P_phase_goal1)).^2);
            
        elseif type==0
            err_eye_freq = ( ( (log(E_gain)-log(E_gain1)) - (log(RR.E_gain_goal0)-log(RR.E_gain_goal1) ) ).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(  (E_phase-E_phase1) - (RR.E_phase_goal0-RR.E_phase_goal1) )).^2 );
            err_pc_freq = ((log(P_gain)- log(RR.P_gain_goal0)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(P_phase - RR.P_phase_goal0)).^2);
            
        elseif type==2
            err_eye_freq = ( ( (log(E_gain)-log(E_gain1)) - (log(RR.E_gain_goal2)-log(RR.E_gain_goal1) ) ).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(  (E_phase-E_phase1) - (RR.E_phase_goal2-RR.E_phase_goal1) )).^2 );
            err_pc_freq = ((log(P_gain)- log(RR.P_gain_goal2)).^2 + ...
                I.fit_phase*deg2rad(wrapTo180(P_phase - RR.P_phase_goal2)).^2);
            
        end
        
        
        % Baseline error: JR data in the dark
        SSE = @(xhat, x) cellfun(@(xhat, x) nansum((xhat - x).^2), xhat, x);
        if type == 1
            mask_JR = ~JR.light & (1:length(JR.light))'<=I.nConds_JR;
            [Ehat_JR, Phat_JR] = getModelJR(KF_X, I, JR.head(mask_JR), JR.target(mask_JR), JR.light(mask_JR), JR.sines(mask_JR));
            
            err_eye_JR = SSE(Ehat_JR, JR.E(mask_JR));
            err_pc_JR = SSE(Phat_JR, JR.P(mask_JR));
        else
            % Specific to learning: only look at Ramachandran & Lisberger frequency data
            err_eye_JR = [];
            err_pc_JR = [];
        end
        
        % STEP
        step_ind = strcmp(JR.conds,'step');
        tt_step = JR.tts{step_ind};
        head_step = JR.head{step_ind}/max(JR.head{step_ind});
        mask_step_ss =  tt_step > I.T_step_ss;  % Start comparing 500 ms in
        mask_step_trans =  tt_step > 0 & tt_step < I.T_step_trans; % Transient is the peak prior to 50 ms - increase to
        
        [Ehat_step, Phat_step] = getModelJR(KF_X, I, {head_step}, {zeros(size(tt_step))}, 0, 0);
        err_eye_step = sum( ( -Ehat_step(mask_step_ss) - WL.eye_ss(type+1)).^2 );
        err_pc_step  = sum( ( Phat_step(mask_step_ss) - WL.pc_ss(type+1)).^2 );
        
        % Get the error in the transient
        if type ~= 1
            temp = Phat_step(mask_step_trans);
            pc_trans_curr = temp(find(abs(temp)==max(abs(temp)),1));
            err_pc_step_trans = (pc_trans_curr - WL.pc_trans(type+1)).^2; % Transient = peak firing in 1st 50 ms after step onset
        else
            % Don't penalize the transient if we are before learning
            err_pc_step_trans  = 0;
        end
        weight_freqs = 100000;
        weight_step = 10000; 
        weight_transient = I.fit_transient*10000; % Don't weight this too much - max measurement of Lisberger 1994 II. may exagerate transient in the data
        
        err_eye = [err_eye_JR;...
            weight_freqs.*err_eye_freq;... 
            weight_step*err_eye_step];              % Step steady state
        err_pc = [err_pc_JR; ...
            weight_freqs.*err_pc_freq; ...
            weight_step*err_pc_step;...
            weight_transient*err_pc_step_trans]; 
        
        
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



end

