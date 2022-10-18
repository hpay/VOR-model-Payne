function [Ehat, Phat, Phat_source] = modelClosedloop(K, I, head_vel, target_vel, ...
    light, sine, pc_stim, N)

% Optional: to simulate a population of PCs, provide a matrix for K.PH and/or K.PR.
% Each col is a cell, rows are filter timepoints
ncells = size(K.PH, 2);

% Optional: add noise in visual pathway
if ~exist('N','var')
    N.noise_PR = 0;
    N.noise_PT = 0;
    noise_on = false;
else
    noise_on = true;
end

% To split up source of PC activity
if nargout >= 3
    return_Psource = true;
else
    return_Psource = false;
end

dt = K.dt;
nt = length(head_vel);
head_vel = head_vel(:);
light = light(:); % light is usually just 0 or 1. It is a vector for occlusion expmt
target_vel = target_vel(:);
target_vel(isnan(target_vel)) = 0;

freqs_JR = I.freqs_JR; % should be [.5 2 5 10]; 

if mean(head_vel==0) > .95 && mean(target_vel==0)>.95
    impulse_test = 1;
else
    impulse_test = 0;
end

% Get the target acceleration exactly right even for sines
if sine && ~impulse_test
    target_acc = diff([target_vel; target_vel(1)])/dt; 
else
    target_acc = diff([target_vel; 0])/dt;
end

% If sine, duplicate vectors to reach steady state
if ~isfield(I,'sine_steadystate'); I.sine_steadystate = 1; end
if sine && ~impulse_test && ~noise_on && I.sine_steadystate % Exception for impulse test
    nrep = round(I.T_steadystate / 2);
    nt = nt*nrep;
    head_vel = repmat(head_vel, nrep, 1);
    target_vel = repmat(target_vel, nrep, 1);
    target_acc = repmat(target_acc, nrep, 1);
end

% Stim to PC should be zero unless specified
if ~exist('pc_stim','var')
    pc_stim = zeros(size(head_vel));
end

% If predictive target input enumerated by freq
if sine 
    T_index = ismember(freqs_JR, sine);
    
    % Interpolate if needed
    if I.include_T_vel
        if any(T_index)
            target_gain_vel = K.PT_vel(T_index,:);
        else
            target_gain_vel = interp1(freqs_JR, K.PT_vel, sine); % TODO switch to log frequency
        end
    end
    
    if  I.include_T_acc
        if any(T_index)
            target_gain_acc = K.PT_acc(T_index,:);
        else
            target_gain_acc = interp1(freqs_JR, K.PT_acc, sine);
        end
    end
    
end

%% Calculate eye and PC timepoint-by-timepoint
[Phat, Ehat, retslip_vel] = deal(zeros(nt,1));
[PR_vel, PT_vel, PT_acc, PR, PT] = deal(zeros(nt,1));
Phat_source = zeros(nt, 3, 2, ncells);
if ncells>1
    [Phat, PR_vel, PR, PT, PT_vel, PT_acc] = deal(zeros(nt, ncells));
    
    % So far, not adding any noise to target weights, but need to keep same
    % size/net weight
    if exist('target_gain_vel','var') && length(target_gain_vel)==1
        target_gain_vel = target_gain_vel*ones(1,ncells)/ncells;
        target_gain_acc = target_gain_acc*ones(1,ncells)/ncells;
    end
end

% Delay for retinal slip (in time stamps)
delayR = find(K.PR_vel, 1);

for ii = 1:nt-1
    
    i_temp_PH = min(ii, size(K.PH,1));
    i_temp_EH = min(ii, size(K.EH,1)); % Make sure this is long enough
    i_temp_EP = min(ii, size(K.EP,1));
    i_temp_PE = min(ii, size(K.PE,1));
    i_temp_vis = min(ii, size(K.PR_vel,1));
    
    % Visual input
    if any(light)
        
        % Calculate proportional and integrated retinal slip
        retslip_vel(ii) = (target_vel(ii) - head_vel(ii)- Ehat(ii));
        
        % Retinal slip velocity
        if length(light) == 1
            PR_vel(ii,:) = retslip_vel(ii:-1:ii-i_temp_vis+1)' *...
                K.PR_vel(1:i_temp_vis,:);
            
            
            
        else
            % Visual occlusion test
            % Only let RS input continue until the RS delay
            PR_vel(ii,:) = retslip_vel(ii:-1:ii-i_temp_vis+1)' * ...
                K.PR_vel(1:i_temp_vis,:);
            
            if ii>delayR && ~light(ii-delayR)
                 PR_vel(ii,:)  = 0;
            end
        end
        
        % ---- Predictive target velocity ----
        PT_vel(ii,:) = 0; % Default

        % Predictive target velocity: sines
        if I.include_T_vel && sine && ii>=(I.delayR/K.dt)            
            if ~impulse_test
                PT_vel(ii,:) = target_vel(ii) * target_gain_vel;
            end            
        end
        
        % Predictive target velocity: steps
        if I.include_T_vel_step && ~sine 
            i_temp_T = min(ii, size(K.PT_vel_step,1));            
            PT_vel(ii,:) = target_vel(ii:-1:ii-i_temp_T+1)' * ...
                K.PT_vel_step(1:i_temp_T,:);
            
        end
        
        % ---- Predictive target acceleration ----
        PT_acc(ii,:) = 0;  % Default
        
        % Predictive target accelleration: sines
        if I.include_T_acc && sine && ii>=(I.delayR/K.dt)
            if  ~impulse_test
                PT_acc(ii,:) = target_acc(ii) * target_gain_acc;
            end
        end
        
        % Predictive target accelleration: steps
        if I.include_T_acc_step && ~sine 
            i_temp_T = min(ii, size(K.PT_acc_step,1));
            PT_acc(ii,:) = target_acc(ii:-1:ii-i_temp_T+1)' * ...
                K.PT_acc_step(1:i_temp_T,:);
            
        end
        
        
        % Net visual input to PC
        PR(ii,:) = PR_vel(ii,:);
        PT(ii,:) = PT_vel(ii,:) + PT_acc(ii,:);
        
        % Add visual noise
        if noise_on
            PR(ii,:) = PR(ii,:)*(1 + N.noise_R(ii,:));
            PT(ii,:) = PT(ii,:)*(1 + N.noise_T(ii,:));
        end
        
    else
        PR(ii,:) = 0;
        PT(ii,:) = 0;
    end
    
    % Head inputs:     
    PH = head_vel(ii:-1:ii-i_temp_PH+1)' * K.PH(1:i_temp_PH, :);
    
    % Eye efference copy:
    if K.PF % For speed, only calc feedback if there is any pos feedback
        PE =  Ehat(ii:-1:ii-i_temp_PE+1)' * K.PE(1:i_temp_PE, :);
    else
        PE = 0;
    end
    
    % Head inputs to eye:
     EH =  head_vel(ii:-1:ii-i_temp_EH+1)' * K.EH(1:i_temp_EH);
    
    % Net PC + visual input to brainstem
%     EP = (K.p * Phat(ii:-1:ii-i_temp_EP+1, :) + ...
%         (1-K.p) * ( PR(ii:-1:ii-i_temp_EP+1, :)+PT(ii:-1:ii-i_temp_EP+1, :) ) )' ...
%         * K.EP(1:i_temp_EP,:);
    EP = Phat(ii:-1:ii-i_temp_EP+1, :)'* K.EP(1:i_temp_EP,:);
    
    % Add vestibular noise if needed
    if noise_on
        PH = PH*(1 + N.noise_PH(ii)) + N.noise_constP(ii); % TODO: change this to just a single constant noise, added at end?
        EH = EH*(1 + N.noise_EH(ii)) + N.noise_constE(ii);
    end
    
    % Purkinje cell update
    Phat(ii+1,:) = PH + PR(ii,:) + PT(ii,:) + PE + pc_stim(ii,:);
    
    % Split up Purkinje cell sources of activity
    if return_Psource

        % Deal with multiple cells!
        Phat_source(ii+1, 1, 1, :) = PH;
        Phat_source(ii+1, 2, 1, :) = PR(ii, :);
        Phat_source(ii+1, 2, 2, :) = PT(ii, :);
        Phat_source(ii+1, 3, 1, :) = PE;
    end
    
    % Eye update
    if ncells == 1        
        Ehat(ii+1) = EH + EP;
    else        
        % If multiple Purkinje cells, combine the their output together by averaging
        Ehat(ii+1) = EH + mean(EP);
    end
    
    if noise_on
        Ehat(ii+1) = Ehat(ii+1)*(1 + N.noise_E(ii));
    end
    
end

% Only report the last cycle of steady state sines
if sine && ~impulse_test && ~noise_on && I.sine_steadystate
    nt = nt/nrep;
    Ehat = Ehat(end-nt+1:end,:);
    Phat = Phat(end-nt+1:end,:);
end
