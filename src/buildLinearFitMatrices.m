function [S, B, S_eye, S_pc] = ...
buildLinearFitMatrices(I, head, target, hevel, PC, sines, light)
% [S, PH_basis, EH_basis, R_basis, T_vel_basis, T_acc_basis,...
% T_vel_step_basis, T_acc_step_basis, PE_basis, EP_basis] = ...
% buildLinearFitMatrices(I, head, target, hevel, PC, sines, light)
%
% Takes the time domain data as input and returns matrices for each
% conditions containing the data multiplied by the basis sets

% Create basis sets
debug_on = 1;
log_scale = 0.1;
norm_on = 1;
B = [];
B.dt = I.dt;
B.B_PH = makeSmoothTemporalBasisCosLog(I.max_time_basis, I.n_basis,I.dt, I.coverage_basis, I.delayH, log_scale, norm_on, debug_on);       % Basis set for PC from head
B.B_EH = makeSmoothTemporalBasisCosLog(I.max_time_basis, I.n_basis,I.dt, I.coverage_basis, I.delayH, log_scale, norm_on, debug_on);       % Basis set for eye from head
B.B_vis = makeSmoothTemporalBasisCosLog(I.max_time_basis_vis, I.n_basis_vis,I.dt, I.coverage_basis, I.delayR, log_scale, norm_on, debug_on);  % Basis set for retinal slip
B.B_EP = makeSmoothTemporalBasisCosLog(I.max_time_basis_EP, I.n_basis,I.dt, I.coverage_basis, I.delayEP, log_scale, norm_on, debug_on);   % Basis set for PC to Eye synapse
B.B_tar = ones(1,length(I.freqs_JR)); % Enumerate - one for each frequency JR tested
B.B_tar_step = 1;

% Create efference copy input
[PE_basis, B.B_PE] = filterExponential(I.tauPE, I.dt, hevel, sines);


% Take target velocity then replace with zeros if dark or sine wave
target_vel_step_zeros = target;
target_vel_step_zeros(~light | sines>0) = cellfun(@(x) zeros(size(x)), target(~light | sines>0),'Uni',0);

% Filter with given tau for predictive target during step
[T_vel_step_basis, B.B_PT_vel_step] = filterExponential(I.T_step_tau_guess, I.dt, target_vel_step_zeros, sines, I.delay_tar_step);

% Do the same thing with target acceleration
target_acc_step_zeros = cellfun(@(x)[0; diff(x)]/I.dt,target_vel_step_zeros,'UniformOutput',false);
[T_acc_step_basis, B.B_PT_acc_step] = filterExponential(I.T_step_tau_guess, I.dt, target_acc_step_zeros, sines, I.delay_tar_step);



% Calculate retinal slip
retslip_vel = cellfun(@(T, H, E) T-H-E, target, head, hevel,'Uni',0);

% Create a predictive target signal that is only present in the light
target_vel = target; % Create copy of target signal, replace dark periods with NaN (for calculating acc and normalizing data)
target_vel(~light) = cellfun(@(x) nan(size(x)), target(~light),'Uni',0);

% Create predictive target acceleration
target_acc = cellfun(@(x)[0; diff(x)]/I.dt,target_vel,'UniformOutput',false);

[PH_basis, EH_basis, R_basis,...
    T_vel_basis, T_acc_basis, EP_basis] = deal(cell(I.nConds,1));

for jj = 1:I.nConds
    
    % Get input history matrices
    pad_data_or_zeros = sines(jj)>0; % Pad sine data with wrapped data, pad steps with zeros
    
    curr_H_vel = makeHistoryMatrix(head{jj}, size(B.B_PH,1), pad_data_or_zeros);
    curr_R_vel = makeHistoryMatrix(retslip_vel{jj}, size(B.B_vis,1), pad_data_or_zeros);
    curr_P = makeHistoryMatrix(PC{jj}, size(B.B_EP,1), pad_data_or_zeros);
    
    % Multiply history matrix with basis set
    PH_basis{jj} = curr_H_vel*B.B_PH;
    EH_basis{jj} = curr_H_vel*B.B_EH;
    R_basis{jj} = curr_R_vel*B.B_vis;
    EP_basis{jj} = curr_P*B.B_EP;

    % Build up target enumeration for each freq
    if sines(jj) && light(jj)
        
        % T_vel_basis{jj}: nt x n_freqs_JR (one for each frequency tested in light)
        % All but one column is zero, correct freq column is filled
        T_vel_basis{jj} = zeros(length(target_vel{jj}), length(I.freqs_JR));
        T_acc_basis{jj} = zeros(length(target_vel{jj}), length(I.freqs_JR));        
        T_vel_basis{jj}(:, sines(jj) == I.freqs_JR)  = target_vel{jj};
        T_acc_basis{jj}(:, sines(jj) == I.freqs_JR)  = target_acc{jj};
        
    else        
        T_vel_basis{jj} = zeros(length(target_vel{jj}), size(B.B_tar,2));
        T_acc_basis{jj} = zeros(length(target_vel{jj}), size(B.B_tar,2));
    end    
    
end


%% Normalize so different inputs have roughly equal amplitude (all are
% already centered around 0). TODO: could eliminate and combine with
% regularization
S = [];
S.scale_P = 1/std(cell2mat(PC(1:I.nConds_JR)),'omitnan');
S.scale_H_vel = 1/std(cell2mat(head(1:I.nConds_JR)));
S.scale_R_vel = 1/std(cell2mat(retslip_vel(1:I.nConds_JR)),'omitnan');
S.scale_T_vel = 1/std(cell2mat(target_vel(1:I.nConds_JR)),'omitnan');
S.scale_T_acc = 1/std(cell2mat(target_acc(1:I.nConds_JR)),'omitnan');
S.scale_E = 1/std(cell2mat(hevel(1:I.nConds_JR)),'omitnan');



%% Define stimulus (S) and response (R) for initial linear regression

% Create stimulus matrix for EYE
S_eye = [cell2mat(EH_basis)*S.scale_H_vel, cell2mat(EP_basis)*S.scale_P];
S_eye(isnan(S_eye)) = 0;

% Define masks for the kernel for each input signal (in basis space)
[B.bEH_maskE, B.bEP_maskE] = deal(false(1,size(S_eye,2)));
order_eye = cumsum([1  I.n_basis]);
B.bEH_maskE(order_eye(1):order_eye(2)-1) = true;
B.bEP_maskE(order_eye(2):end) = true;

% Create stimulus matrix for PC
S_pc = [cell2mat(PH_basis)*S.scale_H_vel, ...                       % Head vel
    cell2mat(R_basis)*S.scale_R_vel ,...                            % Retslip vel
    iif(I.include_T_sine, cell2mat(T_vel_basis)*S.scale_T_vel),...   % Predictable target
    iif(I.include_T_sine, cell2mat(T_acc_basis)*S.scale_T_acc),...
    iif(I.include_T_step, cell2mat(T_vel_step_basis)*S.scale_T_vel),...
    iif(I.include_T_step, cell2mat(T_acc_step_basis)*S.scale_T_acc),...
    cell2mat(PE_basis)...                                           % Efference copy feedback
    ];
S_pc(isnan(S_pc)) = 0;

% Define masks locating the coefficients for each input signal
[B.bPH_maskP, B.bPR_vel_maskP,  B.bPT_vel_maskP,...
    B.bPT_acc_maskP, B.bPE_maskP, B.bPT_vel_step_maskP, B.bPT_acc_step_maskP]  ...
    = deal(false(1,size(S_pc,2)));
order_pc = cumsum([1 I.n_basis,  I.n_basis_vis,...
    I.include_T_sine*size(B.B_tar,2), I.include_T_sine*size(B.B_tar,2),...
    I.include_T_step*[1 1] ]);
B.bPH_maskP(order_pc(1):order_pc(2)-1) = true;
B.bPR_vel_maskP(order_pc(2):order_pc(3)-1) = true;
B.bPT_vel_maskP(order_pc(3):order_pc(4)-1) = true;
B.bPT_acc_maskP(order_pc(4):order_pc(5)-1) = true;
B.bPT_vel_step_maskP(order_pc(5):order_pc(6)-1) = true;
B.bPT_acc_step_maskP(order_pc(6):order_pc(7)-1) = true;
B.bPE_maskP(end) = true;



% Store masks for each filter in a structure
% Mask into "Kb_eye_pc" for all params
B.bEP_mask     = [B.bEP_maskE        false(size(B.bPH_maskP))];
B.bEH_mask     = [B.bEH_maskE        false(size(B.bPH_maskP))];
B.bPH_mask     = [false(size(B.bEH_maskE)) B.bPH_maskP];
B.bPR_vel_mask = [false(size(B.bEH_maskE)) B.bPR_vel_maskP];
B.bPT_vel_mask = [false(size(B.bEH_maskE)) B.bPT_vel_maskP];
B.bPT_vel_step_mask = [false(size(B.bEH_maskE)) B.bPT_vel_step_maskP];
B.bPT_acc_step_mask = [false(size(B.bEH_maskE)) B.bPT_acc_step_maskP];
B.bPT_acc_mask = [false(size(B.bEH_maskE)) B.bPT_acc_maskP];
B.bPE_mask = [false(size(B.bEH_maskE)) B.bPE_maskP];

% Select all weights for FINE TUNING BASELINE in the dark
B.tune_mask_vest_eye_pc = B.bEH_mask |  B.bPH_mask | B.bEP_mask;

% Select all vestibular weights for LEARNING
B.learn_mask_vest_eye_pc = B.bEH_mask |  B.bPH_mask;
if I.fixB % Optional: remove brainstem weights
    B.learn_mask_vest_eye_pc(B.bEH_mask) = 0;
end


