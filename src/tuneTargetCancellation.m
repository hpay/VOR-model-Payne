function [K1_new, K0_new, K2_new] = tuneTargetCancellation(K1, K0, K2, ...
    B, I, S, mask_learn, head, target, conds, light, sines)
% Allow the predictive target signal weights to change
% with learning. Don't change anything about K1 (model before learning)

% Set to 0 to used close loop calculation (for nonlinear retinal slip)
impulse_or_closedloop = 0; 

lb = -inf(nnz(mask_learn),1); % Don't allow negative target signals? CHANGE LB to zero HERE
ub = inf(nnz(mask_learn),1);

% Model helper: impulse
model_helper_impulse = @(F, ind) modelImpulse(F, I, head{ind}, target{ind},...
    light(ind), sines(ind));

% Model helper: slower closed loop version
model_helper_closedloop = @(K, ind) modelClosedloop(K, I, head{ind}, target{ind},...
    light(ind), sines(ind));

K1_new = K1;                    % K1 doesn't change in this model
F1 = getImpulse(K1, I);      % Get the impulse before learning
F1_new = F1;                    % K1 doesn't change in this model

% Original free parameters
p_orig = K1.Kb_eye_pc(mask_learn);


% Fit with fmincon
options = optimoptions('fmincon','MaxIter',100,'FunctionTolerance',0.1,'StepTolerance',1e-5,'MaxFunctionEvaluations', 1e8,...
   'PlotFcn',@optimplotfval, 'Display','iter');
[p2_new, fval, exitflag, output] = fmincon(@(p) myMinFun(p, K2), p_orig,[],[],[],[], lb, ub,[], options);
[p0_new, fval, exitflag, output] = fmincon(@(p) myMinFun(p, K0), p_orig,[],[],[],[], lb, ub,[], options);


K2_new = updateKp(K2, p2_new, mask_learn);
K0_new = updateKp(K0, p0_new, mask_learn);


    function [err_out, eye_canc_normal, eye_canc_learn] = myMinFun(p, K)
        
        % Update strucutre K with new variable weights
        Kx = updateKp(K, p, mask_learn);
        
        if impulse_or_closedloop
            % Get impulse response
            Fx = getImpulse(Kx, I);
            
            % 1b. Calculate error during cancellation, compared to baseline
            [err_eye_canc, eye_canc_normal, eye_canc_learn] = get_canc_error_impulse(F1, Fx);
        else
            
            % 1b. Calculate error during cancellation, compared to baseline
            [err_eye_canc, eye_canc_normal, eye_canc_learn] = get_canc_error_closedloop(K1, Kx);
        end
        
        % Avoid the baseline drifting
        err_baseline = mean(eye_canc_learn).^2; % TODO  check compared to mean(eye_canc_learn.^2)
     
        % Original error
        err_out =  err_eye_canc+err_baseline;
   
    end

    function [err_canc_learn, eye_canc_normal, eye_canc_learn] = get_canc_error_impulse(F1, FX, ploton)
        if ~exist('ploton','var')
            ploton = 0;
        end
        
        % Cancellation index
        ind_canc = strcmp(conds,'05Hz_x0');
        
        % Eye during cancellation before learning
        eye_canc_normal = model_helper_impulse(F1, ind_canc);
        
        % Eye during cancellation after "learning":
        eye_canc_learn = model_helper_impulse(FX, ind_canc) ;
        
        % Error for difference in the eye during cancellation
        err_canc_learn = mean((eye_canc_normal - eye_canc_learn).^2);
        
        if ploton
            cla;
            plot(eye_canc_normal,'k'); hold on; plot(eye_canc_learn,'m');
            legend('Normal', 'Learn'); title('Eye cancellation')
        end
    end


    function [err_canc_learn, eye_canc_normal, eye_canc_learn] = get_canc_error_closedloop(K1, KX, ploton)
        if ~exist('ploton','var')
            ploton = 0;
        end
        
        % Cancellation index
        ind_canc = strcmp(conds,'05Hz_x0');
        
        % Eye during cancellation before learning
        eye_canc_normal = model_helper_closedloop(K1, ind_canc) ;
        
        % Eye during cancellation after "learning":
        eye_canc_learn = model_helper_closedloop(KX, ind_canc) ;
        
        % Error for difference in the eye during cancellation
        err_canc_learn = mean((eye_canc_normal - eye_canc_learn).^2);
        
        if ploton
            cla;
            plot(eye_canc_normal,'k'); hold on; plot(eye_canc_learn,'m');
            legend('Normal', 'Learn'); title('Eye cancellation')
        end
    end


    function K_out = updateKp(K_in, p_new, mask_learn)
        Kb_eye_pc_new = K_in.Kb_eye_pc;
        Kb_eye_pc_new(mask_learn) = p_new;
        
        
        Kb_eye = Kb_eye_pc_new(1:length(K_in.Kb_eye));
        Kb_pc = Kb_eye_pc_new(length(K_in.Kb_eye)+1:end);
        
        K_out = updateK(B, I, S,  Kb_eye, Kb_pc);
    end

end