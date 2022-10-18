
function [Ehat, Phat] = modelImpulse(F, I, head_vel, target_vel, light, sine, debug_on)

head_vel = head_vel(:);
target_vel = target_vel(:);


if ~exist('debug_on','var')
    debug_on= 0;
end

freqs_JR = [.5 2 5 10];

nt = length(head_vel);
 
% Create target acc and pos in case needed for predictive target signal
target_vel(isnan(target_vel)) = 0;

target_acc = diff([target_vel; target_vel(1);])/F.dt;

% Duplicate vectors to reach steady state
if sine && I.sine_steadystate
    
    if isfield(I,'T_steadystate')
        nrep = round(I.T_steadystate / 2);
    else
        nrep = 3;
    end
    
    head_vel = repmat(head_vel, nrep, 1);
    target_vel = repmat(target_vel, nrep, 1);
    target_acc = repmat(target_acc, nrep, 1);
      
end

% Calculate eye and PC
if ~light
    % HEAD INPUT - IN THE DARK
    E_head = filter(F.EH, 1, head_vel);
    P_head = filter(F.PH, 1, head_vel);
    
    % VISUAL INPUT - IN THE DARK
    E_vis = zeros(size(head_vel));
    P_vis = zeros(size(head_vel));
    
else

    if sine
            
        % HEAD INPUT - IN THE LIGHT - different from dark!
%         E_head = filter(F.EH_light_sine, 1, head_vel);
%         P_head = filter(F.PH_light_sine, 1, head_vel);       
        E_head = filter(F.EH_light, 1, head_vel);
        P_head = filter(F.PH_light, 1, head_vel);        
        
        % VISUAL INPUT - IN THE LIGHT
        E_vis = filter(F.ET_sine, 1, target_vel);
        P_vis = filter(F.PT_sine, 1, target_vel);
        
        % If the target option = 'enum'
        mask = freqs_JR==sine;
        if length(F.PT_vel)==length(mask)
            
            % Separate step for retinal slip visual input and predictive
            % target input
            P_target_vel = F.PT_vel(mask)*target_vel;
            P_target_acc = zeros(size(target_vel));
            
            if length(F.PT_acc)==length(mask); P_target_acc = F.PT_acc(mask)*target_acc;
            end
 
            P_target = P_target_vel + P_target_acc;
            
            % Filter by the impulse response of a stim to PCs (due to network
            % feedback)
            E_sine = filter(F.EP_light, 1, P_target);
            P_sine = filter(F.PP_light, 1, P_target);
            
            E_vis = E_vis + E_sine;
            P_vis = P_vis + P_sine;
            
        end
        
        
    else % step
        
                    
        % HEAD INPUT - IN THE LIGHT - different from dark!
%         E_head = filter(F.EH_light_step, 1, head_vel);
%         P_head = filter(F.PH_light_step, 1, head_vel);   
        
        E_head = filter(F.EH_light, 1, head_vel);
        P_head = filter(F.PH_light, 1, head_vel);        
        
        % VISUAL INPUT - IN THE LIGHT
        E_vis = filter(F.ET_step, 1, target_vel);
        P_vis = filter(F.PT_step, 1, target_vel);
        
        
        
    end
    
    
end

% Net output: Y_total = Y_head + Y_vis
Ehat = E_head + E_vis;
Phat = P_head + P_vis;

% DEBUG
if debug_on
    figure;
    subplot(211); plot(E_head); hold on; plot(E_vis); title('Eye: head and vis')
    subplot(212); plot(P_head); hold on;  plot(P_vis); title('PC: head and vis')
end

% Only include the second cycle if sine
if sine && I.sine_steadystate
    Ehat = Ehat(end-nt+1:end,:);
    Phat = Phat(end-nt+1:end,:);
end




