function  fH = plotInstability(I,  K, c)


% Generate input signals
T_smooth = 0.025;   % (s) % 0.01
tbefore = -0.1;  
tafter = 0.25;      % (s) Time to go out for step
amp = 10;           % (deg/s) Scale for step

dt = K.dt;
tt = tbefore:dt:tafter;

% Generate smoothed step
head_curr = double(tt>0)*amp;
head_curr = smooth(head_curr, round(T_smooth/dt));
head_curr = smooth(head_curr, round(T_smooth/dt));
% K, I, head_vel, target_vel,  light, sine, pc_stim, N)


% Change weights
K_PH = K;
K_EH = K;
K_PH.PH = K_PH.PH*1.1;
K_EH.EH = K_EH.EH*1.1;

% Run model
[Ehat_orig, Phat_orig]     = modelClosedloop(K, I, head_curr, zeros(size(head_curr)), 0, 0);
[Ehat_PH, Phat_PH] = modelClosedloop(K_PH, I, head_curr, zeros(size(head_curr)), 0, 0);
[Ehat_EH, Phat_EH] = modelClosedloop(K_EH, I, head_curr, zeros(size(head_curr)), 0, 0);


% Plot result
plot(tt, Ehat_orig,'Color',c); hold on
plot(tt, Ehat_PH,'--', 'Color',c); hold on
plot(tt, Ehat_EH,':', 'Color',c); hold on

