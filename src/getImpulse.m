%% 7a. Find impulse for closed loop model
function F = getImpulse(K, I, option_dark_only)
% 
if ~exist('option_dark_only','var')
    option_dark_only = 0;
end

F = [];
F.dt = K.dt;

% Predictive target signal includes position -- can't calculate impulse
% response correctly (?), so just pass on the info for now
F.PT_vel = K.PT_vel;
F.PT_acc = K.PT_acc;

% if option_dark_only
%     T = 2; % (s)
% else
    T = I.T_steadystate;
% end
nt = round(T/K.dt); % Due to feedback loops needs to be as long as the longest simulation

%==========1. HEAD: in dark ========
% Define length of net filter to be twice single filter (since output from
% PC will then get filtered to become output of eye)
delta = zeros(nt,1); delta(1) = 1; % Impulse
null = zeros(nt,1);

% Simulate params - head input only
sine_on = 0;
light_on = 0;
head_curr = delta; % Impulse
target_curr = null;

% Run model
[Ehat, Phat] = modelClosedloop(K, I, head_curr, target_curr, ...
    light_on, sine_on); % light on, sine on
F.EH = Ehat(1:end);         % EH impulse
F.PH = Phat(1:end);         % PH impulse

if option_dark_only
    
    F.EH_light = 0;
    F.PH_light = 0;
    
    F.ET_step = 0;
    F.PT_step = 0;
    
    F.ET_sine = 0;
    F.PT_sine = 0;
    
    F.PP_light = 0;
    F.EP_light = 0;        
    
else
    
    %========== HEAD with lights on, sine  ========
    sine_on = 1;
    light_on = 1;
    head_curr = delta; % Impulse
    target_curr = null;
    
    [Ehat_light, Phat_light] = modelClosedloop(K, I, head_curr, target_curr, ...
        light_on, sine_on); % light on, sine on
    F.EH_light = Ehat_light(1:end);         % EH impulse
    F.PH_light = Phat_light(1:end);         % PH impulse
    
%     %========== HEAD with lights on, step  ========
%     sine_on = 0;
%     light_on = 1;
%     head_curr = delta; % Impulse
%     target_curr = null;
%     
%     [Ehat_light, Phat_light] = modelClosedloop(K, I, head_curr, target_curr, ...
%          light_on, sine_on); % light on, sine on
%     F.EH_light_step = Ehat_light(1:end);         % EH impulse
%     F.PH_light_step = Phat_light(1:end);         % PH impulse
    
    
    % ============= VISUAL: sine===============
    sine_on = 1;
    light_on = 1;
    head_curr = null;
    target_curr = delta;
    [Ehat_light_sine, Phat_light_sine] = modelClosedloop(K, I, head_curr, target_curr, ...
        light_on, sine_on); % light on, sine on
    F.ET_sine = Ehat_light_sine(1:end);
    F.PT_sine = Phat_light_sine(1:end);
    
    
    % ===== VISUAL: step============
    sine_on = 0;
    light_on = 1;
    head_curr = null;
    target_curr = delta;
    [Ehat_light_step, Phat_light_step] = modelClosedloop(K, I, head_curr, target_curr, ...
        light_on, sine_on); % light on, sine on
    F.ET_step = Ehat_light_step(1:end);
    F.PT_step = Phat_light_step(1:end);
    
        
    %=========5. Input to PC, in the light ============    
    light_on = 1;
    sine_on = 0;
    head_curr = null;
    target_curr = null;
    PC_stim_curr = delta;
    
    [Ehat, Phat] = modelClosedloop(K, I, head_curr, target_curr, ...
        light_on, sine_on, PC_stim_curr); % light on, sine on
    F.EP_light = Ehat(1:end);
    F.PP_light = Phat(1:end);
    
end



% Eliminate long ends of some filters with zeros
 filter_mask = @(x)  (1:length(x)) <= find( diff([1; abs(x); 0]<1e-6 )>0, 1, 'last');
 F.EH = F.EH(filter_mask(F.EH));
 F.PH = F.PH(filter_mask(F.PH));
 
 
 F.EH_light = F.EH_light(filter_mask(F.EH_light));
 F.PH_light = F.PH_light(filter_mask(F.PH_light));
%  F.EH_light_sine = F.EH_light_sine(filter_mask(F.EH_light_sine));
%  F.PH_light_sine = F.PH_light_sine(filter_mask(F.PH_light_sine));
%  F.EH_light_step = F.EH_light_step(filter_mask(F.EH_light_step));
%  F.PH_light_step = F.PH_light_step(filter_mask(F.PH_light_step));
 F.ET_step = F.ET_step(filter_mask(F.ET_step));
 F.PT_step = F.PT_step(filter_mask(F.PT_step));
 F.ET_sine = F.ET_sine(filter_mask(F.ET_sine));
 F.PT_sine = F.PT_sine(filter_mask(F.PT_sine));
 F.EP_light = F.EP_light(filter_mask(F.EP_light));
 F.PP_light = F.PP_light(filter_mask(F.PP_light));

end