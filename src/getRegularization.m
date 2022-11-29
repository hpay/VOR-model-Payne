function R = getRegularization(B, I)
% ************ Eye fit: first to get scaling for kPE (kEP*kPE = PF)**********
R = struct;
R.reg_weights_eye = [I.scale_EH*(1:nnz(B.bEH_maskE)).^I.reg_exp ...
    I.scale_EP*(1:nnz(B.bEP_maskE)).^I.reg_exp_EP];

% Get mask where coeffients switch from one filter to another
n_eye = length(B.bEH_maskE);
mask_edge_eye = false(n_eye,1);
mask_edge_eye([find(B.bEH_maskE,1) find(B.bEP_maskE,1) ]) = true;

% Get Tikhonov regularization matrices
R.T0_eye = makeTikhonov(R.reg_weights_eye, 0);
R.T1_eye = makeTikhonov(R.reg_weights_eye, 1, mask_edge_eye);
R.T2_eye = makeTikhonov(R.reg_weights_eye, 2, mask_edge_eye);


% Sum of weights, weight increasing with time
R.reg_weights_pc = [I.scale_PH*(1:nnz(B.bPH_maskP)).^I.reg_exp ...
    I.scale_PR*(1:nnz(B.bPR_vel_maskP)).^I.reg_exp...
    I.scale_PT*[ones(1, nnz(B.bPT_vel_maskP))...
    ones(1, nnz(B.bPT_acc_maskP))]...
    iif(I.include_T_step, I.scale_PT*[1 1])...
    0 ];

% Get mask where coefficients switch from one filter to another
n_pc = length(B.bPH_maskP);
mask_pc_edge = false(n_pc,1); % Don't penalize curvature between edges of different filters
mask_pc_edge([find(B.bPH_maskP,1)  find(B.bPR_vel_maskP,1)...
    find(B.bPT_vel_maskP)  find(B.bPT_acc_maskP) ...
    find(B.bPT_vel_step_maskP) find(B.bPT_acc_step_maskP) find(B.bPE_maskP)   ]) = true;

% Get Tikhonov regularization matrices
R.T0_pc = makeTikhonov(R.reg_weights_pc, 0);
R.T1_pc = makeTikhonov(R.reg_weights_pc, 1, mask_pc_edge);
R.T2_pc = makeTikhonov(R.reg_weights_pc, 2, mask_pc_edge);

