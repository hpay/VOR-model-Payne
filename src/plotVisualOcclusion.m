function plotVisualOcclusion(I,  K, c)

% Parameters
amp = 10;           % (deg/s)
t_before = 0.2;     % (s)
t_after = 1;        % (s)
t_light_off = 0.5;  % (s)
t_light_on = 1;     % (s)
t_smooth = 0.02;

% Run simulation
tt = -t_before:K.dt:t_after; 
head_curr = zeros(size(tt));
target_curr = smooth(amp*(tt>0), round(t_smooth/K.dt));
target_curr = smooth(target_curr, round(t_smooth/K.dt));
light = ones(size(tt));
light(tt>t_light_off & tt<t_light_on) = 0;
sine = 0;

[Ehat, Phat] = modelClosedloop(K, I, head_curr, target_curr, ...
    light, sine);

% Simulate the same thing but without the blink
% [Ehat_normal, Phat_normal] = modelClosedloop(K, I, head_curr, target_curr, ...
%     ones(size(tt)), sine);

% Plot results
ax1 = gca;
% Plot target and eye 
h = plot(tt, target_curr + amp,'k'); hold on
plot(tt(tt>t_light_off), target_curr(tt>t_light_off) + amp,'k:','Clipping','off')
h(end+1) = plot(tt, Phat+amp/2,'Color', c,'LineWidth',1); hold on
h(end+1) = plot(tt, Ehat,'Color', c,'LineWidth',2); hold on
% h(end+1) = plot(tt, Ehat_normal+amp/2,'Color', c,'LineWidth',2); hold on
% h(end+1) = plot(tt, Ehat,'Color', c,'LineWidth',2); hold on

fixticks
xlabel('Time (s)')
ax1.YAxis.Visible = 'off'; % remove y-axis
ax1.XAxis.Visible = 'off'; % remove y-axis
set(ax1,'XLim', [-t_before t_after+.2],'YLim',[-2 23])

% Plot scale bars
x = 1.05; y = 2;
dx = 0.2; dy = 10;
plot(x+[0 dx], y*[1 1], 'k','Clipping','off')
plot((x+dx)*[1 1], y+[0 dy],'k','Clipping','off')

% Plot a shaded box
ha = rectangle('Pos',[t_light_off min(ylim) t_light_on-t_light_off range(ylim)], 'FaceColor',.8*[1 1 1],'EdgeColor','none');
uistack(ha,'bottom')
uistack(h(1),'bottom')

% legend(h, {'Target','PC','Eye'},'box','off','Position',[.7 .6 .3 .2])


