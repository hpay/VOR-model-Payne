function [conds, tts, head, target, hevel, PC, sines, light] = addStep(I, conds, ...
    tts, head, target, hevel, PC, sines, light)

% Add an ideal VORD step to enforce steady state
tt_step = (1:round(I.T_step/I.dt))'*I.dt + I.T_step_start; % Step time vector
step_in = I.scale_step * double(tt_step>0);
step_in = filter(ones(round(I.tau_step/I.dt),1),round(I.tau_step/I.dt), step_in);  % Smooth it
step_in = filter(ones(round(I.tau_step/I.dt),1),round(I.tau_step/I.dt), step_in);  % Smooth it
e_step = I.eye_gain_ss*step_in;     % Define desired eye steady state
p_step = I.pc_gain_ss*step_in;      % Define desired PC steady state

conds = [conds; 'step'];            % Tack on this ideal step input to all the relevant variables
tts = [tts; {tt_step}];
head = [head; {step_in}];
target = [target; {zeros(size(step_in))}];
hevel = [hevel; {e_step}];
PC = [PC; {p_step}];
sines = [sines; 0];
light = [light; 0];


