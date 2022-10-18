function [K0new, K1new, K2new] = tuneVestibularSimultaneous(...
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

% Combine all Jennifer's data in one structure
JL.conds = conds;
JL.tts = tts;
JL.head = head;
JL.target = target;
JL.E = E;
JL.P = P;
JL.light = light;
JL.sines = sines;

% Goal gain and phase of eye
RL.E_gain_goal1 = RL_data.gains1/RL_data.scaleJR2SL;
RL.E_phase_goal1 = RL_data.phases1;
RL.E_gain_goal2 = RL_data.gains2/RL_data.scaleJR2SL;
RL.E_phase_goal2 = RL_data.phases2;
RL.E_gain_goal0 = RL_data.gains0/RL_data.scaleJR2SL;
RL.E_phase_goal0 = RL_data.phases0;

% Get the goal PC baseline gain and phase for 0.5 Hz (close to steady state)
[pc_amp, pc_phase] = fitsine(dt, P{1}, .5);

vord_step_ind = strcmp(conds,'500ms_dark');
pc_gain_ss = I.pc_gain_ss;
% pc_gain_ss = (nanmean(P{vord_step_ind}(400:600))/nanmean(head{vord_step_ind}(400:600)));
trans_mask = tts{vord_step_ind}>0 & tts{vord_step_ind}< I.T_step_trans;
temp1 = P{vord_step_ind}(trans_mask);
pc_trans_init = temp1(find(abs(temp1)==max(abs(temp1)),1))/max(head{vord_step_ind});

WL.eye_ss = [ RL.E_gain_goal0(1) RL.E_gain_goal1(1) RL.E_gain_goal2(1)];

% Combining SL and Watanabe
WL.pc_ss = pc_gain_ss + 1/RL_data.scaleJR2SL*...
    I.pc_ss.*(1-I.eye_ss); % change in sp/s per change in deg/s, * change in gain


% Create a PC goal gain and phase
RL.P_gain_goal1 = NaN(size(RL_data.freqs));
RL.P_phase_goal1 = NaN(size(RL_data.freqs));
RL.P_gain_goal1(1) = abs(pc_gain_ss);
RL.P_phase_goal1(1) = pc_phase;

RL.P_gain_goal0 = NaN(size(RL_data.freqs));
RL.P_phase_goal0 = NaN(size(RL_data.freqs));
RL.P_gain_goal0(1) = abs(WL.pc_ss(1));
RL.P_phase_goal0(1) = 11; % deg, Watanabe 1985

RL.P_gain_goal2 = NaN(size(RL_data.freqs));
RL.P_phase_goal2 = NaN(size(RL_data.freqs));
RL.P_gain_goal2(1) = abs(WL.pc_ss(3));
RL.P_phase_goal2(1) = 172; % deg, Watanabe 1985

RL.freqs = RL_data.freqs;

% Using Lisberger data for the transients (don't weight too much, likely to
% ohigh because it measures the peak rate in 1st 50 ms)
WL.pc_trans =  pc_trans_init + 1/RL_data.scaleJR2SL*[1.62  0 -1.01];% TODO: check the scaling

% Original free parameters (3 x n)
p_init = [K0_0.Kb_eye_pc(mask_tune)...
    K1_0.Kb_eye_pc(mask_tune)...
    K2_0.Kb_eye_pc(mask_tune)];

% Find the initial net weight of kEP to enforce that it doesn't change
% sumEP = sum(K1_0.EP); % EDIT 8/22/19
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

[p_out, fval, exitflag, output] = fmincon(@(x)myMinFun(x, B, I, S,...
    K1_0, K2_0, K0_0, mask_tune, mask_fixEP, JL, RL, WL, R), p_init,A_ineq,B_ineq,...
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
p1(mask_fixEP) = p1(mask_fixEP) * sumEP_orig/sumEP_final;
p0(mask_fixEP) = p1(mask_fixEP);
p2(mask_fixEP) = p1(mask_fixEP);

K0new = updateKp(K1_0, B, I, S, p0, mask_tune);
K1new = updateKp(K1_0, B, I, S, p1, mask_tune);
K2new = updateKp(K1_0, B, I, S, p2, mask_tune);

    function err_out = myMinFun(p, B, I, S, K1_0, K2_0, K0_0, mask_tune, mask_fixEP, JL, RL, WL, R)
        tic
        p0_curr = p(:,1);
        p1_curr = p(:,2);
        p2_curr = p(:,3);
        
        % Enfore kEP isn't changing before learning: Force total weight of kEP to stay the same
        sumEP = sum(K1_0.EP);
        K1_0_init = updateKp(K1_0, B, I, S,  p1_curr, mask_tune);        % Temporary variable
        sumEP_curr = sum(K1_0_init.EP)
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
        [err_eye_all0, err_pc_all0, err_reg0] = getErrorVest(I, p0_curr, p1_curr,...
            KF0,KF1,0, mask_tune, JL, RL, WL, R);
        [err_eye_all1, err_pc_all1, err_reg1] = getErrorVest(I, p1_curr, p1_curr,...
            KF1,KF1,1, mask_tune, JL, RL, WL, R);
        [err_eye_all2, err_pc_all2, err_reg2] = getErrorVest(I, p2_curr, p1_curr,...
            KF2,KF1,2, mask_tune, JL, RL, WL, R);
        
        % Sum error for each learning condition
        err_out0 = [nansum(err_eye_all0)  nansum(err_pc_all0) I.scale_learn*err_reg0];
        err_out1 = [nansum(err_eye_all1)  nansum(err_pc_all1)  err_reg1];
        err_out2 = [nansum(err_eye_all2)  nansum(err_pc_all2)  I.scale_learn*err_reg2];
        
        % Error for deviation from original EP weight % TODO: check necessary?
        err_EP  = I.scale_const_EP * (sumEP_curr - sumEP).^2;
        
        % Total error
        err_out_curr = [err_out0  err_out1  err_out2  err_EP];
%         legend_str = {'E0','P0','reg0','E1','P1','reg1','E2','P2','reg2','scaleEP'};
        err_out = sum(err_out_curr);
        toc
    end




end

