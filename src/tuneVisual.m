function [K1_new, K0_new, K2_new] = tuneVisual(K1, K0, K2, ...
    B, I, S, mask_learn, head, target, E, P,  conds, light, sines, R)
% Tune the visual weights (retinal slip and target)  for all models at the same
% time. Fit the performance of eye and PC on all visual conditions (17 conds)
% before learning, and optionally minimize change in cancellation performance over learning
%
% I.option_weight_05Hz_x0 = 1; % Set to weight for even weighting of JR's data, change to 10 to emphasize x0 condition at 0.5 Hz
% I.option_weight_500ms = 1; % Set to weight for even weighting of JR's data, change to 10 to emphasize x0 condition at 0.5 Hz
%
%  error = sum_conds( (y-yhat)^2 ) + beta*(y_can_x0 - y_canc_x1)^2 + beta*(y_canc_x2 - y_canc_x1)^2;

close all
dt = K1.dt;

if I.impulse_or_closedloop
    % Model helper: impulse
    modelHelper = @(F, ind) modelImpulse(F, I, head{ind}, target{ind},...
        light(ind), sines(ind));
else
    % Model helper: slower closed loop version
    modelHelper = @(K, ind) modelClosedloop(K, I, head{ind}, target{ind},...
        light(ind), sines(ind));
end
% Original free parameters
p_orig = K1.Kb_eye_pc(mask_learn);


lb = -inf(length(K1.Kb_eye_pc),1);
ub = inf(length(K1.Kb_eye_pc),1);

% If enumerating predictive target, the tau for steps (last param) must be positive:
if I.fit_T_step_tau~=0
    lb(find(B.bPT_vel_step_mask, 1, 'last')) = 0;
    lb(find(B.bPT_acc_step_mask, 1, 'last')) = 0;
end

% If all target inputs must be positive
if I.restrict_PT_pos
    lb(B.bPT_vel_step_mask) = 0;
    lb(B.bPT_acc_step_mask) = 0;
    lb(B.bPT_vel_mask) = 0;
    lb(B.bPT_acc_mask) = 0;
end

lb = lb(mask_learn);
ub = ub(mask_learn);


% METHOD 1: fminsearch (slower,  but better for non-smooth function)
options = optimset('PlotFcns',@optimplotfval,'TolFun',10,'TolX',1e-4, 'MaxIter',1e8); % Tol fun from 0.1
[p_new, fval, exitflag, output] = fminsearchbnd(@(p) myMinFun(p), p_orig, lb, ub, options);

% METHOD 3: fmincon (faster for smooth functions)
% options = optimoptions('fmincon','MaxIter',100,'FunctionTolerance',1,...
%     'StepTolerance',1e-6,'MaxFunctionEvaluations', 1e8,...
%    'PlotFcn',@optimplotfval, 'Display','iter','Algorithm','sqp');
% [p_new, fval, exitflag, output] = fmincon(@myMinFun, p_orig,[],[],[],[], lb, ub,[], options);

disp(output.message)

K0_new = updateKp(K0, p_new, mask_learn);
K1_new = updateKp(K1, p_new, mask_learn);
K2_new = updateKp(K2, p_new, mask_learn);

    function err_out = myMinFun(p)
        n_weights = nnz(mask_learn);
        
        % Update strucutre K with new variable weights
        K1_curr = updateKp(K1, p, mask_learn);
        K2_curr = updateKp(K2, p, mask_learn);
        K0_curr = updateKp(K0, p, mask_learn);
        
        if I.impulse_or_closedloop
            % Get impulse response
            KF1_curr = getImpulse(K1_curr, I);
        else
            KF1_curr = K1_curr;
        end
        
        
        % Calculate error in baseline visual conds for eye & PC
        [err_eye_all, err_pc_all, Ehat, Phat] = getVisualError(KF1_curr);       
        
        % Sum eye and pc error
        err_eye_vis = sum(err_eye_all);
        err_pc_vis = sum(err_pc_all);
        
        % Create combined regularization matrices for both eye and PC params
        Tall0  = blkdiag(R.T0_eye, R.T0_pc);
        Tall1  = blkdiag(R.T1_eye, R.T1_pc);
        Tall2  = blkdiag(R.T2_eye, R.T2_pc);
        
        % Mask out only the params we care about now
        Tall0_mask = Tall0(mask_learn, mask_learn);
        Tall1_mask = Tall1(mask_learn, mask_learn);
        Tall2_mask = Tall2(mask_learn, mask_learn);
        
        % Regularization penalty:
        err_reg = sum(I.reg_lambda0*(Tall0_mask*p(1:n_weights)).^2 +...
            I.reg_lambda1*(Tall1_mask*p(1:n_weights)).^2 +...
            I.reg_lambda2*(Tall2_mask*p(1:n_weights)).^2);
        
        % Penalty for drift in steady state
        
        % Total error
        err_out = err_eye_vis + err_pc_vis + err_reg;
        
        if isnan(err_out)
            pause(.1)
        end
        
        % Split error up by source
        % err_out_split(1,:) = [err_eye_vis err_pc_vis  scale_canc_error*err_eye_canc_low scale_canc_error*err_eye_canc_high];
    end

    function [err_canc_learn, eye_canc_normal, eye_canc_learn] = getCancError(KF1, KFx, ploton)
        if ~exist('ploton','var')
            ploton = 0;
        end
        
        % Cancellation index
        ind_canc = strcmp(conds,'05Hz_x0');
        
        % Eye during cancellation before learning
        eye_canc_normal = modelHelper(KF1, ind_canc);
        
        % Eye during cancellation after "learning":
        eye_canc_learn = modelHelper(KFx, ind_canc) ;
        
        % Error for difference in the eye during cancellation
        err_canc_learn = sum((eye_canc_normal - eye_canc_learn).^2);
        
        if ploton
            cla;
            plot(eye_canc_normal,'k'); hold on; plot(eye_canc_learn,'m');
            legend('Normal', 'Learn'); title('Eye cancellation')
        end
    end

% Get error before learning for eye and PC, for all visual conditions
    function [err_eye, err_pc, E_light_hat, P_light_hat] = getVisualError(KF)
        light_inds = find(light);
        E_light = E(light_inds);
        P_light = P(light_inds);
        [err_eye, err_pc] = deal(NaN(1,length(light_inds)));
        [P_light_hat, E_light_hat] = deal(cell(1, length(light_inds)));
        parfor ii_temp = 1:length(light_inds)
            curr_ind = light_inds(ii_temp);
            [E_light_hat{ii_temp}, P_light_hat{ii_temp}] = modelHelper(KF, curr_ind);
            err_eye(ii_temp) = sum((E_light_hat{ii_temp} - E_light{ii_temp}).^2);
            err_pc(ii_temp) = sum((P_light_hat{ii_temp} - P_light{ii_temp}).^2);
        end
    end


    function K_out = updateKp(K_in, p_new, mask_learn)
        Kb_eye_pc_new = K_in.Kb_eye_pc;
        
        Kb_eye_pc_new(mask_learn) = p_new(1:nnz(mask_learn));
        
        Kb_eye = Kb_eye_pc_new(1:length(K_in.Kb_eye));
        Kb_pc = Kb_eye_pc_new(length(K_in.Kb_eye)+1:end);
        
        K_out = updateK(B, I, S, Kb_eye, Kb_pc);
        
        
    end

end