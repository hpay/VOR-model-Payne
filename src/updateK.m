function K_out = updateK(B, I, S,  Kb_eye, Kb_pc)
% Update all the relevant things in the K structure given a new set of eye
% and PC weights (concatenated)

K_out = [];
K_out.dt = B.dt;
K_out.p = I.p; 

Kb_eye_pc_new = [Kb_eye(:); Kb_pc(:)];

% TODO: fix this assignment if implementing saturation
% K_out.PR_vel_sat = I.saturateR_vel;

    
% Default values from optional inputs
K_out.EH = [];
K_out.EP = [];

K_out.PH = [];

K_out.PR_vel = [];

K_out.PT_vel = 0;
K_out.PT_acc = 0;
K_out.PT_vel_step = 0;
K_out.PT_acc_step = 0;

%%%% EYE %%%%
    
% To Eye from Head:
K_out.EH = S.scale_H_vel * B.B_EH * Kb_eye_pc_new(B.bEH_mask);

% To Eye from PC:
K_out.EP = S.scale_P * B.B_EP * Kb_eye_pc_new(B.bEP_mask);

%%%% PC %%%%

% To PC from Eye:
K_out.PE = B.B_PE * Kb_eye_pc_new(B.bPE_mask);

% To PC from Head:
K_out.PH = S.scale_H_vel * B.B_PH * Kb_eye_pc_new(B.bPH_mask);

% To PC from visual:
K_out.PR_vel = S.scale_R_vel * B.B_vis * Kb_eye_pc_new(B.bPR_vel_mask);

% Predictive target input
K_out.PT_vel = S.scale_T_vel * Kb_eye_pc_new(B.bPT_vel_mask);
K_out.PT_acc = S.scale_T_acc * Kb_eye_pc_new(B.bPT_acc_mask);

if I.include_T_vel_step
    
    % Create exponential filter with given delay
    PT_vel_step_gain = Kb_eye_pc_new(find(B.bPT_vel_step_mask,1));
    
    if I.fit_T_step_tau
        PT_vel_step_tau = Kb_eye_pc_new(find(B.bPT_vel_step_mask,1,'last'));
    else
        PT_vel_step_tau = I.T_step_tau_guess;
    end
    
    [~, B_PE] = filterExponential(PT_vel_step_tau, B.dt, {1}, 0, I.delay_tar_step);
    K_out.PT_vel_step = B_PE;
    
    % Scale according to scale factor and gain
    K_out.PT_vel_step = S.scale_T_vel * PT_vel_step_gain * K_out.PT_vel_step;
    
end


if I.include_T_acc_step
    
    % Create exponential filter with given delay
    PT_acc_step_gain = Kb_eye_pc_new(find(B.bPT_acc_step_mask,1));
    
    if I.fit_T_step_tau
        PT_acc_step_tau = Kb_eye_pc_new(find(B.bPT_acc_step_mask,1,'last'));
    else
        PT_acc_step_tau = I.T_step_tau_guess;
    end
        
    [~, B_PE] = filterExponential(PT_acc_step_tau, B.dt, {1}, 0, I.delay_tar_step);
    K_out.PT_acc_step = B_PE;
    
    % Scale according to scale factor and gain
    K_out.PT_acc_step =  S.scale_T_acc * PT_acc_step_gain * K_out.PT_acc_step;
    
end
    

K_out.Kb_eye_pc = Kb_eye_pc_new;
K_out.Kb_eye = Kb_eye;
K_out.Kb_pc = Kb_pc;

% Check pos feedback (round to nearest thousandth)
K_out.PF = round(sum(K_out.PE) * sum(K_out.EP)*1e3)/1e3;
      

end