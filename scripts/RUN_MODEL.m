function RUN_MODEL(savename)
% RUN_MODEL Fit of closed loop model parameters to visual and vestibular data
%   RUN_MODEL(filepathname) runs all feedback strengths: 0, 0.1, ... 1 and
%   stores the results in the folder filepathname within results
%   RUN_MODEL(FILEPATHNAME, FEEDBACKS) runs the specified strengths
%   (FEEDBACK = sum(kEP)*sum(kPE))
%
%   For example, run the No-Feedback and Strong-Feedback models:
%       RUN_MODEL('C:\Dropbox\rlab\model\TO_SHARE\20190821', [0 1])
%
%   Dependencies: All files in the subfolder helperFunctions, and the subfolder
%   data
%
%   Author: Hannah Payne, Goldman & Raymond labs, 2019
%
%
% =============================DATASETS====================================
%   Dataset 1: Jennifer Raymond's recordings of PURKINJE CELL and EYE VELOCITY
%   in response to visual and vestibular inputs -- steps and
%   sines. Specified in the model as conds, tts, head, target, hevel, PC
%
%   Dataset 2: Ramachandran & Lisberger 2005 EYE VELOCITY data in response to
%   vestibular inputs in the dark, from 0.5 to 50 Hz, in Normal monkeys and
%   after High and Low gain training. Specified in the model as RR_data.
%
%   Dataset 3: Steady state/low frequency PURKINJE CELL responses to VOR
%   before v. after learning (Watanabe 1985, Lisberger 1994b). Specified in
%   the model as I.eye_ss: eye steady state after learning and I.pc_ss:
%   Purkinje cell steady state after learning
%
%
% =============================INSTRUCTIONS=============================
% Run RUN_MODEL.m from within a folder named "VOR-model".
%
% If RUN_MODEL.m is located in e.g. ~\model_and_results\VOR-model\RUN_MODEL.m,
% results will be placed in the next folder up, within a subfolder created
% with current date and iteration: e.g. ~\model_and_results\20190801_1
% Existing results folders will not be overwritten.
%
% To replot existing results, run PLOT_MODEL.m with the correct folder
% path as the argument .
%
% See comments below for descriptions of parameters
%
% =============================RESULTS=============================
% For each positive feedback strength specified
%
% I: structure of parameters
% S: structure of scale factors (stdev or raw signal) used to normalize
% data before fitting
% B: basis set
% R:
% RR_data: structure containing data quantified from Ramachandran and
% Lisberger 2005, their Figure
% K1: baseline fit before learning
% K2: fit after VOR-increase learning
% K0: fit after VOR-decrease learning
% conds: list of behavioral conditions used to fit model
% E: error of fits
%
% After all positive feedback strengths have been run, a final results file
% will be saved containing all the individual fits combined (for variables K1, K0, K2, and E).


% Add helper function folders to search path
code_folder = fileparts(fileparts(which(mfilename)));
addpath(genpath(fullfile(code_folder, 'src')))

warning off

%% =================== PARAMETERS ==================

I.PFs = (0:.1:1)';       % Positive feedback strengths. Generally don't change, but can change to negative to test neg feedback

% Specify the index of feedback strengths to run here. They will be run in
% this order.
I.runPFs = I.PFs;

% Parameters - data loading
I.data_path = fullfile(code_folder, 'data');       % Where is data located?
I.t_step_jr = .7;                               % (s) The max time to display after start of step
I.RR_data_mean_or_W = 'meanfix';                % {'W','mean','Wscale','meanfix'} which monkey from Ramachandran Lisberger 2005 to pick %*** 9/11
I.data_add_interp_PC = true;                    % Interpolate JR data to predict PC responses between 0.5 and 10 Hz for frequency fine tuning

I.option_save_results = 1;                      % Save pdfs
I.option_save_PFs = [-1 -.5 -.9 0 0.2 0.5 0.8 0.9 1];   % Save traces from which models?
I.impulse_or_closedloop = 1;    % Set 1 to use impulse response: faster, but everything must be linear

% Parameters - fitting, predictive target signal
I.include_T_vel = 1;            % {0,1} Include predictive target velocity for sines
I.include_T_acc = 1;            % {0,1} Include predictive target acceleration for sines
I.include_T_vel_step = 1;       % {0,1} Include predictive target velocity for steps
I.include_T_acc_step = 1;       % {0,1} Include predictive target acceleration for steps
I.T_step_tau_guess = 0.025;      % (s) Guess for tau of acc and vel for step prediction. This will get fine tuned if next option is set.
I.fit_T_step_tau = 0;           % {0,1} Set 0 to not fit the above tau [***TODO: currently not working. But 50 ms is a good guess***]
I.tunePT = 1;                   % {0,1} allow predictive target signal weight to be different for x0, x1, and x2 learning conditions

% Parameters - fitting, predictive target signal
% I.include_T_vel = 0;            % {0,1} Include predictive target velocity for sines
% I.include_T_acc = 0;            % {0,1} Include predictive target acceleration for sines
% I.include_T_vel_step = 0;       % {0,1} Include predictive target velocity for steps
% I.include_T_acc_step = 0;       % {0,1} Include predictive target acceleration for steps
% I.T_step_tau_guess = 0.025;      % (s) Guess for tau of acc and vel for step prediction. This will get fine tuned if next option is set.
% I.fit_T_step_tau = 0;           % {0,1} Set 0 to not fit the above tau [***TODO: currently not working. But 50 ms is a good guess***]
% I.tunePT = 0;                   % {0,1} allow predictive target signal weight to be different for x0, x1, and x2 learning conditions

% Paremeters - fitting, general
I.eye_ss = [.4; 1; 1.6];     % VOR gain to achieve after learning
I.pc_ss = [ (1.36 + .72)/2;  0; (.75 + .55)/2]; % change in sp/s per change in deg/s reported by Lisberger 1994II and Watanabe 1985, averaged together
I.sine_steadystate = 1;         % {0,1} Take steady state sine response (keep last cycle only)
I.T_steadystate = 6;            % (s) Time at which to measure "steady state". Must be a multiple of 2 s!
I.scale_const_EP = 1e6;         % (big number) Scale factor to enforce that the EP synapse only changes shape, not net weight, when fine tuning
I.amp = 10;                     % (deg/s) scale artificial sine waves to match stim amplitude in JR's data
I.p = 1; % [0 1] (default 1) Fraction of retinal slip to go through PC pathway -- rest is rerouted. NOT CURRENTLY IMPLEMENTED
I.algorithm = 'interior-point'; % Default for linear fits

% Include long artificial step input to ensure steady state responses are correct
I.include_step = 1;             % {0,1} Include an idealized step in order to fit steady state eye and PC before and after learning
I.scale_step = 100;             % Weight the step more -- changed from 100 to 1000 1/5/2019. Changed back 8/20
I.T_step = 2;                   % (s) Duration of step
I.T_step_ss = 0.5;              % (s) Measure steady state after this time
I.T_step_trans = 0.1;           % (s) Measure max transient before this time
I.T_step_start = -0.1;          % (s) Time before step input begins
I.tau_step = 0.025;             % (s) Smooth step

% Parameters for basis set - include string "basis" somewhere!
I.max_time_basis = .05;         % (s) How long out to go in time: PH, EH
I.max_time_basis_EP = .15;      % (s) How long out to go in time: EP
I.max_time_basis_vis = 0.5;     % (s) How long out to go in time: retinal slip pathway  
I.n_basis = 10;                 % Number of bases % ***12
I.n_basis_vis = 12;             % Number of bases for retinal slip
I.coverage_basis = 2;           % Set >1 to allow overlap of adjacent basis functions (default 2)
I.fit_phase = 1; % TODO IMPLEMENT these
I.fit_transient = 1; 
%  TODO: add back?   I.option_weight_500ms = 10; % Set to weight for even weighting of JR's data, change to 10 to emphasize 500 ms cond
% I.scale_vis_error = 1;  % How much to weight visual error
%     I.scale_canc_error = 0; % How much to weight cancellation error: eye should be similar before v. after learning
%
% Parameters: time constants and delays
I.delayH = 0.005;               % (s) Min delay for head input filter
I.delayR = 0.060;               % (s) Min delay for retinal slip input filter % changed from .05
I.delayEP = 0.001;              % (s) Min delay for Purkinje cell->brainstem synapse filter
I.delay_tar_step = I.delayR;    % (s) Min delay for visual target filter (steps only)
I.tauPE = .003;                 % (s) Time constant of efference feedback (fixed)

% Parameters: Constrain weights to test necessity of certain things
I.fixB = 0;                     % Change to 1 to block any brainstem plasticity
I.fixPH = 0;                    % Change to 1 to only allow LTP (no LTD), -1 only allow LTD (no LTP)
I.fixPH_absolute = 0;           % TODO: IMPLEMENT *** Change to 1 to only allow LTP (no LTD), -1 only allow LTD (no LTP) (ANY WEIGHT)
I.restrict_PT_pos = 1;          % constrain predictive target weights positive
I.restrict_EP_pos = 0;          % constrain PC-> brainstem synapse to be all positive

% Parameters: regularization
I.reg_exp = 2;                  % Exponent for increasing regularization penalty with coefficient number
I.reg_exp_EP = 0;               % Same as above for EP synapse

I.scale_PH = 1;          % Scale regularization penalty for PH pathway (default 1)
I.scale_EH = 1;          % Scale regularization penalty for EH pathway (default 1)
I.scale_EP = 40;         % Scale regularization penalty for EP pathway
I.scale_PR = 6;          % relative scale for retinal slip input 
I.scale_PT = 0.25;       % relative scale for predictive target inputs 

I.reg_lambda0 = 1;       % scale regularization penalty for sum squared weight of model coefficients
I.reg_lambda1 = 0;       % scale derivative penalty (slope of basis set coefficients) (default 0)
I.reg_lambda2 = 1;       % scale second derivative penalty (curvature of basis set coefficients) (default 0)
I.scale_learn = 12;      % Scale the regularization penalty for changes in kPH and kEH during learning

% Colors
I.c = [linspace(46,247,length(I.PFs))' linspace(49, 148, length(I.PFs))' linspace(146, 30, length(I.PFs))']/255;  % PF strength
I.learn_color = [.7*[1 1 1];.4*[1 1 1];[0 0 0]]; % Learning condition

I.savename = savename;
pathname = fullfile(code_folder, 'results');
I.figures_path = fullfile(pathname, I.savename);
disp(I.figures_path);

% Create a copy of the current code to keep track of changes
mkdir(I.figures_path)
copyfile(fullfile(code_folder,'src'), fullfile(I.figures_path,'src'))
copyfile(fullfile(code_folder,'scripts'), fullfile(I.figures_path,'scripts'))

%% =================== LOAD INPUTS ===================
[conds, I.nConds_JR, tts, head, target, hevel, PC, sines, light,...
    I.dt, n_cells, RR_data] = loadJR_RR_combined(I, I.data_path);

% Measure "steady state" eye and PC  from JR data
vord_step_ind = strcmp(conds,'500ms_dark');
ss_mask = tts{vord_step_ind}>=0.15 & tts{vord_step_ind}<0.25;
I.pc_gain_ss = (nanmean(PC{vord_step_ind}(ss_mask))/nanmean(head{vord_step_ind}(ss_mask)));
I.eye_gain_ss = (nanmean(hevel{vord_step_ind}(ss_mask))/nanmean(head{vord_step_ind}(ss_mask)));

% Add an ideal VORD step to enforce steady state
[conds, tts, head, target, hevel, PC, sines, light] = addStep(I, conds, ...
    tts, head, target, hevel, PC, sines, light);
I.nConds = length(conds);

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
I.freqs_JR = unique(sines(1:I.nConds_JR)); I.freqs_JR(I.freqs_JR==0) = [];
B.B_tar = ones(1,length(I.freqs_JR)); % Enumerate - one for each frequency JR tested
if  ~I.fit_T_step_tau;   B.B_tar_step = 1;
else B.B_tar_step = [1 1]; end


% Create stimulus history matrices in basis space
[S, B, PH_basis, EH_basis, R_basis, T_vel_basis, T_acc_basis,...
    T_vel_step_basis, T_acc_step_basis, PE_basis, EP_basis] = ...
    buildLinearFitMatrices(I, B, head, target, hevel, PC, sines, light);

disp('Bases multiplied')

%% 1. Define stimulus (S) and response (R) for initial linear regression

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
    iif(I.include_T_vel, cell2mat(T_vel_basis)*S.scale_T_vel),...   % Predictable target
    iif(I.include_T_acc, cell2mat(T_acc_basis)*S.scale_T_acc),...
    iif(I.include_T_vel_step, cell2mat(T_vel_step_basis)*S.scale_T_vel),...
    iif(I.include_T_acc_step, cell2mat(T_acc_step_basis)*S.scale_T_acc),...
    cell2mat(PE_basis)...                                           % Efference copy feedback
    ];
S_pc(isnan(S_pc)) = 0;

% Define masks locating the coefficients for each input signal
[B.bPH_maskP, B.bPR_vel_maskP,  B.bPT_vel_maskP,...
    B.bPT_acc_maskP, B.bPE_maskP, B.bPT_vel_step_maskP, B.bPT_acc_step_maskP]  ...
    = deal(false(1,size(S_pc,2)));
order_pc = cumsum([1 I.n_basis,  I.n_basis_vis,...
    I.include_T_vel*size(B.B_tar,2), I.include_T_acc*size(B.B_tar,2),...
    I.include_T_vel_step*size(B.B_tar_step,2), I.include_T_acc_step*size(B.B_tar_step,2)]);
B.bPH_maskP(order_pc(1):order_pc(2)-1) = true;
B.bPR_vel_maskP(order_pc(2):order_pc(3)-1) = true;
B.bPT_vel_maskP(order_pc(3):order_pc(4)-1) = true;
B.bPT_acc_maskP(order_pc(4):order_pc(5)-1) = true;
B.bPT_vel_step_maskP(order_pc(5):order_pc(6)-1) = true;
B.bPT_acc_step_maskP(order_pc(6):order_pc(7)-1) = true;
B.bPE_maskP(end) = true;

% Define desired response (R) for linear regression
R_pc = cell2mat(PC);        % Desired pc
R_eye = cell2mat(hevel);    % Desired eye

% Only look at low frequencies (which has both eye and PC data) for now (fine tune later)
sine_mask = repelem(sines, cellfun(@length,head));
lo_freq_mask = sine_mask<=10; % TODO is this necessary?
S_pc_lo = S_pc(lo_freq_mask,:);
S_eye_lo = S_eye(lo_freq_mask,:);
R_pc_lo = R_pc(lo_freq_mask);
R_eye_lo = R_eye(lo_freq_mask);

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

%% Initialization before running loop
tic

%% 2. ===================MAIN LOOP THROUGH POSITIVE FEEDBACK VALUES===================
for ii_temp = 1:length(I.runPFs)
    
    rng(0)

    ii = find(abs(I.PFs-I.runPFs(ii_temp))<eps);
    close all;
    fprintf('Starting PF = %g\n', I.PFs(ii));
        
            
    % Avoid re-running
    curr_savename = sprintf('%s_%gPF.mat', I.savename, I.PFs(ii));    
    curr_savepathname = fullfile(I.figures_path, curr_savename);
    if exist(curr_savepathname,'file')
        continue;
    end
    %% 1. LINEAR REGRESSION TO INITIALIZE WEIGHTS
    
    % ************ Eye fit: first to get scaling for kPE (kEP*kPE = PF)**********
    R = struct;
    R.reg_weights_eye = [I.scale_EH*(1:nnz(B.bEH_maskE)).^I.reg_exp ...
        I.scale_EP*(1:nnz(B.bEP_maskE)).^I.reg_exp_EP];
    
    % Get mask where coeffients switch from one filter to another
    mask_edge_eye = false(size(S_eye_lo,2),1);
    mask_edge_eye([find(B.bEH_maskE,1) find(B.bEP_maskE,1) ]) = true;
    
    % Get Tikhonov regularization matrices
    R.T0_eye = makeTikhonov(R.reg_weights_eye, 0);
    R.T1_eye = makeTikhonov(R.reg_weights_eye, 1, mask_edge_eye);
    R.T2_eye = makeTikhonov(R.reg_weights_eye, 2, mask_edge_eye);
    
    % Add regularization penalty
    X_eye = [S_eye_lo; I.reg_lambda0*R.T0_eye; I.reg_lambda1*R.T1_eye; I.reg_lambda2*R.T2_eye]; % Penalize curvature
    Y_eye = [R_eye_lo; zeros(3*size(S_eye_lo,2),1)];
    
    % Constrain head velocity filters
    beq_eye = zeros(size(S_eye_lo,2), 1);
    Aeq_eye = zeros(size(S_eye_lo, 2),1);
    
    % Upper and lower bounds
    lb_eye = -inf(size(S_eye_lo,2), 1);
    ub_eye = inf(size(S_eye_lo,2), 1);
    
    if I.restrict_EP_pos
        lb_eye(B.bEP_maskE) = 0;
    end
    
    % Solve for eye coefficients
    Kb_eye = lsqlin(X_eye, Y_eye,[],[], diag(Aeq_eye), beq_eye, lb_eye, ub_eye, [], optimset('display','off','Algorithm','interior-point'));
    
    % g1 is the net strength of K_EP. Use to scale K_PE so that K_EP*K_PE = I.PFs(ii).
    K_EP = B.B_EP*Kb_eye(B.bEP_maskE)*S.scale_P;
    g1 = sum(K_EP); % g1*g2 = I.PFs(ii) --> g2 = I.PFs(ii)/g1
    
    %  ******* Purkinje cell fit *******
    % Sum of weights, weight increasing with time
    R.reg_weights_pc = [I.scale_PH*(1:nnz(B.bPH_maskP)).^I.reg_exp ... 
        I.scale_PR*(1:nnz(B.bPR_vel_maskP)).^I.reg_exp...
        I.scale_PT*[ones(1, nnz(B.bPT_vel_maskP))...
        ones(1, nnz(B.bPT_acc_maskP))]...
        iif(I.include_T_vel_step, I.scale_PT*[1  iif(I.fit_T_step_tau, 0)])... % First param for kPT_vel_step is gain, 2nd is tau
        iif(I.include_T_acc_step, I.scale_PT*[1  iif(I.fit_T_step_tau, 0)])... % First param for kPT_vel_step is gain, 2nd is tau
        0 ];
    
    % Get mask where coefficients switch from one filter to another
    mask_pc_edge = false(size(S_pc_lo,2),1); % Don't penalize curvature between edges of different filters
    mask_pc_edge([find(B.bPH_maskP,1)  find(B.bPR_vel_maskP,1)...
        find(B.bPT_vel_maskP)  find(B.bPT_acc_maskP) ...
        find(B.bPT_vel_step_maskP) find(B.bPT_acc_step_maskP) find(B.bPE_maskP)   ]) = true;
    
    % Get Tikhonov regularization matrices
    R.T0_pc = makeTikhonov(R.reg_weights_pc, 0);
    R.T1_pc = makeTikhonov(R.reg_weights_pc, 1, mask_pc_edge);
    R.T2_pc = makeTikhonov(R.reg_weights_pc, 2, mask_pc_edge);
    
    X_reg = [I.reg_lambda0*R.T0_pc; I.reg_lambda1*R.T1_pc;  I.reg_lambda2*R.T2_pc]; % EDIT 6/25 remove I.scale_PH
    Y_reg = zeros(3*size(S_pc_lo,2),1);
    
    X_pc = [S_pc_lo; X_reg];
    Y_pc = [R_pc_lo; Y_reg];
    
    % Fix the positive feedback weight to I.PFs(ii)
    Aeq_pc = double(B.bPE_maskP);
    beq_pc = zeros(size(S_pc_lo,2), 1);
    beq_pc(B.bPE_maskP) = I.PFs(ii)/g1; % g1 is the sum of weights kEP
    
    % Set the predictive target time constant for steps. 
    if I.fit_T_step_tau~=0
        Aeq_pc(find(B.bPT_vel_step_maskP, 1, 'last')) = 1; % Fix the tau only
        beq_pc(find(B.bPT_vel_step_maskP, 1, 'last')) = I.T_step_tau_guess;
        Aeq_pc(find(B.bPT_acc_step_maskP, 1, 'last')) = 1; % Fix the tau only
        beq_pc(find(B.bPT_acc_step_maskP, 1, 'last')) = I.T_step_tau_guess;
    end
    
    % Lower and upper bounds
    lb_pc = -inf(size(S_pc_lo,2),1);
    ub_pc = inf(size(S_pc_lo,2),1);
    
    if I.restrict_PT_pos  % constrain predictive target weights positive
        lb_pc(B.bPT_vel_maskP | B.bPT_acc_maskP | B.bPT_vel_step_maskP | B.bPT_acc_step_maskP ) = 0;
    end
    
    % For PF = 1 only, don't allow predictive target velocity for steps (causes instability)
    if I.PFs(ii) == 1 && I.include_T_vel_step
        lb_pc(B.bPT_vel_step_maskP) = 0;
        ub_pc(B.bPT_vel_step_maskP) = 0;
    end
    
    if ii_temp == 1
        
        % Solve for PC coefficients
        Kb_pc = lsqlin(X_pc, Y_pc,[],[], diag(Aeq_pc), beq_pc, lb_pc, ub_pc, [], optimset('display','off','Algorithm', I.algorithm));
        
        fprintf('Linear fit complete for PF = %g\n', I.PFs(ii));
        
        % Multiply coefficients by bases to get filters in time domain
        K1_orig = updateK(B, I, S, Kb_eye, Kb_pc);
        
        % Plot linear fit predictions of eye and PC
        [err_eye, err_pc, Ehat_linear, Phat_linear] = plotLinearPredictions(S_eye_lo, S_pc_lo,R_eye_lo, R_pc_lo,  Kb_eye, Kb_pc);
        fprintf('Baseline linear fit RMSE, eye: %.3g deg/s,  PC: %.3g sp/s\n', err_eye, err_pc)
        
        % Plot eye frequency response
        freqs = logspace(log10(.5), log10(50), 100);
        [~, ~, E_gain, E_phase, P_gain, P_phase] = getFreq(K1_orig, I, freqs);
        hF_freq = figure;
        [h_freq, hs_freq] = plotFreq(freqs, E_gain*RR_data.scaleJR2SL, E_phase, '-o',I.learn_color(2,:),I.learn_color(2,:));
        subplot(211); title('Eye')
        
        % Plot resulting filter
        figure;  [hf_filter, ~] = plotFilters(K1_orig);
        
        % Plot closed loop model
        E = struct;
        [hf_basefit, E.rmse_eye,  E.rmse_pc, E.nrmse_eye, E.nrmse_pc] = ...
            plotBaselineResults(K1_orig, I, conds, tts, head, target, hevel, PC, light, sines);
        
        if I.option_save_results && ismember(I.PFs(ii), I.option_save_PFs)
            my_export_fig(hf_basefit, fullfile(I.figures_path, sprintf('baseFit%g_orig.pdf', I.PFs(ii))));
            my_export_fig(hf_filter, fullfile(I.figures_path, sprintf('baseFilters%g_orig.pdf', I.PFs(ii))));
        end
        
        % Test models: impulse response and closed loop should match.
        %     TODO delete
            %{
        I.impulse_or_closedloop = 0;
        [hf_basefit, E.rmse_eye,  E.rmse_pc, E.nrmse_eye, E.nrmse_pc] = ...
            plotBaselineResults(K1_orig, I, conds, tts, head, target, hevel, PC, light, sines, [1 2 3 4 23 24 25]);
        
        I.impulse_or_closedloop = 1;
        [hf_basefit, E.rmse_eye,  E.rmse_pc, E.nrmse_eye, E.nrmse_pc] = ...
            plotBaselineResults(K1_orig, I, conds, tts, head, target, hevel, PC, light, sines, [1 2 3 4 23 24 25]);
        %}
        
        % SEQUENTIAL LINEAR FIT - intial guess for learned weights
        K2_orig = tuneVestibularLinear(K1_orig, B, I, S, B.learn_mask_vest_eye_pc, 2, conds, tts,...
            R, RR_data, S_pc, S_eye, R_pc, R_eye);
        
        K0_orig = tuneVestibularLinear(K1_orig, B, I, S, B.learn_mask_vest_eye_pc, 0, conds,tts, ...
            R, RR_data, S_pc, S_eye, R_pc, R_eye);
        
        %  Plot after learning
        [hF_freq_eye, hF_freq_pc, hF_filters, hF_step, hF_sine] = ...
            plotPostLearning(I, K1_orig, K2_orig, K0_orig, RR_data, tts, head, conds);
        
    else
        prev_savename = sprintf('%s_%gPF.mat', I.savename, I.PFs(ii)-0.1);

        P = load(fullfile(I.figures_path, prev_savename));
        K1_orig = P.K1;
        K0_orig = P.K0;
        K2_orig = P.K2;
        
        % Update the feedback loop strength
        g1 = sum(K1_orig.EP); % g1*g2 = I.PFs(ii) --> g2 = I.PFs(ii)/g1
        K1_orig.Kb_pc(B.bPE_maskP) = I.PFs(ii)/g1;        
        K0_orig.Kb_pc(B.bPE_maskP) = I.PFs(ii)/g1;           
        K2_orig.Kb_pc(B.bPE_maskP) = I.PFs(ii)/g1;
        
        K1_orig = updateK(B, I, S, K1_orig.Kb_eye, K1_orig.Kb_pc);
        K0_orig = updateK(B, I, S, K0_orig.Kb_eye, K0_orig.Kb_pc);        
        K2_orig = updateK(B, I, S, K2_orig.Kb_eye, K2_orig.Kb_pc);
        
        % For PF = 1 only, don't allow predictive target velocity for steps (causes instability)
        if I.PFs(ii) == 1 && I.include_T_vel_step
            K1_orig.Kb_pc(B.bPT_vel_step_maskP) = 0;           
            K2_orig.Kb_pc(B.bPT_vel_step_maskP) = 0;            
            K0_orig.Kb_pc(B.bPT_vel_step_maskP) = 0;                       
        end
        
        K1_orig = updateK(B, I, S, K1_orig.Kb_eye, K1_orig.Kb_pc);
        K0_orig = updateK(B, I, S, K0_orig.Kb_eye, K0_orig.Kb_pc);        
        K2_orig = updateK(B, I, S, K2_orig.Kb_eye, K2_orig.Kb_pc);
        
    
    end
    
    
    
    %% 2. NONLINEAR OPTIMIZATION FOR VESTIBULAR WEIGHTS ON CLOSED LOOP MODEL
    % Fine tune vestibular weights and implement vestibular learning
    close all;
    % SIMULTANEOUS FINE TUNING AND LEARNING OF VESTIBULAR PARAMETERS
    fprintf('Starting vestibular fine tuning, simultaneous for PF = %g\n', I.PFs(ii));
    [K0, K1_vest, K2] =  tuneVestibularSimultaneousOld(...
        K0_orig, K1_orig, K2_orig, B, I, S,...
        tts,  head, target, hevel, PC,  conds, light, sines, R, RR_data);
    fprintf('Finished vestibular fine tuning, simultaneous for PF = %g\n', I.PFs(ii));
    
    %     tic
    %     % SIMULTANEOUS FINE TUNING AND LEARNING OF VESTIBULAR PARAMETERS
    %     fprintf('Starting vestibular fine tuning, simultaneous for PF = %g\n', I.PFs(ii));
    %     [K0, K1_vest, K2] =  tuneVestibularSimultaneous(...
    %         K0_temp, K1_orig, K2_temp, B, I, S,...
    %         tts,  head, target, hevel, PC,  conds, light, sines, R, RR_data);
    %     fprintf('Finished vestibular fine tuning, simultaneous for PF = %g\n', I.PFs(ii));
    %     toc
    
    % Plot closed loop model
    [hf_basefit, E.rmse_eye,  E.rmse_pc, E.nrmse_eye, E.nrmse_pc]= ...
        plotBaselineResults(K1_vest, I, conds(~light), tts(~light), head(~light), target(~light), hevel(~light), PC(~light), light(~light), sines(~light));
    
    % Plot for learning
    [hF_freq_eye, hF_freq_pc, hF_filters, hF_step, hF_sine] = ...
        plotPostLearning(I, K1_vest, K2, K0, RR_data, tts, head, conds);
    
    % Plot freq data
    hF_freq_data = figure;
    plotFreq(RR_data.freqs, RR_data.gains1, RR_data.phases1+180, 'o-', I.learn_color(2,:), I.learn_color(2,:));
    plotFreq(RR_data.freqs, RR_data.gains2, RR_data.phases2+180, 'o-', I.learn_color(3,:), I.learn_color(3,:));
    plotFreq(RR_data.freqs, RR_data.gains0, RR_data.phases0+180, 'o-', I.learn_color(1,:), I.learn_color(1,:));
    xlabel('Frequency'); subplot(211); title('Eye')
    shrink([.4 .7]);  fixticks; drawnow
    
    % Save figures
    if I.option_save_results && ismember(I.PFs(ii), I.option_save_PFs)
        my_export_fig(hF_freq_eye, fullfile(I.figures_path, sprintf('learn_freq%g.pdf', I.PFs(ii))));
        my_export_fig(hF_filters, fullfile(I.figures_path, sprintf('learn_filters%g.pdf', I.PFs(ii))),'-dpdf');
        my_export_fig(hF_step, fullfile(I.figures_path, sprintf('learn_step%g.pdf', I.PFs(ii))),'-dpdf');
        my_export_fig(hF_freq_data, fullfile(I.figures_path, sprintf('learn_freq_data.pdf')),'-dpdf');
    end
    
    % Update K1 with tuned weights
    K1 = K1_vest;
    
    
    %% 3. NONLINEAR OPTIMIZATION FOR VISUAL WEIGHTS ON CLOSED LOOP MODEL
    % Fine tune retinal slip and target parameters
    
    
    I.impulse_or_closedloop = 0; % Temp, actually faster this way
    
    % Select the weights to fine tune
    if I.PFs(ii) == 1 && I.include_T_vel_step
        
        % Don't allow predictive target velocity for steps for PF = 1, leads to ramp
        learn_mask_vis_eye_pc = [false(size(B.bEH_maskE)) ...
            (B.bPR_vel_maskP  ...
            | B.bPT_vel_maskP | B.bPT_acc_maskP ...
            | B.bPT_acc_step_maskP)];
    else
        learn_mask_vis_eye_pc = [false(size(B.bEH_maskE)) ...
            (B.bPR_vel_maskP ...
            | B.bPT_vel_maskP | B.bPT_acc_maskP  ...
            | B.bPT_vel_step_maskP  | B.bPT_acc_step_maskP)];
    end
    % Fine tune retinal slip and predictive target weights (not learning, applies to all)
    fprintf('Starting retinal slip fine tuning for PF = %g\n', I.PFs(ii));
    [K1_vis, K0_vis, K2_vis] = ...
        tuneVisual(K1, K0, K2, B, I, S, learn_mask_vis_eye_pc,...
        head, target, hevel, PC,  conds, light, sines, R);
    fprintf('Finished retinal slip fine tuning for PF = %g\n', I.PFs(ii));
    
    
    %  Fine tune changes in predictive target parameters for 0.5 Hz cancellation
    if I.tunePT && (I.include_T_vel || I.include_T_acc)
        
        % Select the weights to allow to change
        learn_mask_vis_eye_pc = [false(size(B.bEH_maskE)) false(size(B.bPH_maskP))];
        learn_mask_vis_eye_pc(length(B.bEH_maskE) + [find(B.bPT_vel_maskP,1) find(B.bPT_acc_maskP,1)]) = true;
        
        %  Tune predictive target signals
        [K1_vis, K0_vis, K2_vis] = ...
            tuneTargetCancellation(K1_vis, K0_vis, K2_vis, B, I, S, learn_mask_vis_eye_pc, head, target, conds, light, sines);
        fprintf('Finished predictive target fine tuning for PF = %g\n', I.PFs(ii));
    end
    
    I.impulse_or_closedloop = 1;
    
    
    %% Visual fine tuning plots
    % Plot new model before learning: with tuned retinal slip/target
    [hf_basefit_vis, E.rmse_eye_vis, E.rmse_pc_vis, ...
        E.nrmse_eye_vis, E.nrmse_pc_vis] =...
        plotBaselineResults(K1_vis, I, conds, tts, head, target, hevel, PC, light, sines);
    
    % 4d. Plot resulting visual filters before v. after learning side by side
    hf_filters = plotFilters(K1_vis);
    
    % Test cancellation for each learning case -- before and after
    hf_canc = plotCancellationVisTuning(...
        I,  K1, K2, K0, K1_vis, K2_vis, K0_vis);
    
    % Save results for this PF strength
    if I.option_save_results && ismember(I.PFs(ii), I.option_save_PFs)
        my_export_fig(hf_basefit_vis, fullfile(I.figures_path, sprintf('baseFit%g_tuneR.pdf', I.PFs(ii))),'-dpdf');
        my_export_fig(hf_filters, fullfile(I.figures_path, sprintf('baseFilters%g_tuneR.pdf', I.PFs(ii))),'-dpdf');
        my_export_fig(hf_canc, fullfile(I.figures_path, sprintf('cancellation_PF%g_tuneR.pdf', I.PFs(ii))),'-dpdf');
    end
    
    %% 8. Save results
    K1 = K1_vis;
    K0 = K0_vis;
    K2 = K2_vis;
    
    save(curr_savepathname, 'I', 'S', 'B', 'R', 'RR_data',...
        'K1','K2','K0','conds',  'E')
    
    fprintf('Complete PF = \n');
    disp(I.PFs(ii))
end

%% Combine all files from this run
combineFiles(I.figures_path)

disp('Files combined and saved')

%% PLOT RESULTS
PLOT_MODEL(I.figures_path)

