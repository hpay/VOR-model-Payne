function [S, B, PH_basis, EH_basis, R_basis, T_vel_basis, T_acc_basis,...
T_vel_step_basis, T_acc_step_basis, PE_basis, EP_basis] = ...
buildLinearFitMatrices(I, B, head, target, hevel, PC, sines, light)
% [S, PH_basis, EH_basis, R_basis, T_vel_basis, T_acc_basis,...
% T_vel_step_basis, T_acc_step_basis, PE_basis, EP_basis] = ...
% buildLinearFitMatrices(I, head, target, hevel, PC, sines, light)
%
% Takes the time domain data as input and returns matrices for each
% conditions containing the data multiplied by the basis sets


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

% Normalize so different inputs have roughly equal amplitude (all are
% already centered around 0). TODO: could eliminate and combine with
% regularization
S = [];
S.scale_P = 1/std(cell2mat(PC(1:I.nConds_JR)),'omitnan');
S.scale_H_vel = 1/std(cell2mat(head(1:I.nConds_JR)));
S.scale_R_vel = 1/std(cell2mat(retslip_vel(1:I.nConds_JR)),'omitnan');
S.scale_T_vel = 1/std(cell2mat(target_vel(1:I.nConds_JR)),'omitnan');
S.scale_T_acc = 1/std(cell2mat(target_acc(1:I.nConds_JR)),'omitnan');
S.scale_E = 1/std(cell2mat(hevel(1:I.nConds_JR)),'omitnan');
