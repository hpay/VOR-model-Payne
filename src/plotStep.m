function plotStep(K_1, K_2, K_0, I, learn_color)
%% plotStep(K_1, K_2, K_0, I, learn_color)
% 
% Plots step response before and after learning

T_smooth = 0.02;   % (s) % 0.01
T_step = 0.25;      % (s) Time to go out for step
amp = 10;           % (deg/s) Scale for step

dt = K_1.dt;
tt_step = (-.05:dt:T_step)';

% Generate smoothed step
head_curr = double(tt_step>0)*amp;
head_curr = smooth(head_curr, round(T_smooth/dt));
head_curr = smooth(head_curr, round(T_smooth/dt));

[Ehat_orig, Phat_orig]     = modelClosedloop(K_1, I, head_curr, zeros(size(head_curr)), 0, 0);
[Ehat_learn2, Phat_learn2] = modelClosedloop(K_2, I, head_curr, zeros(size(head_curr)), 0, 0);
[Ehat_learn0, Phat_learn0] = modelClosedloop(K_0, I, head_curr, zeros(size(head_curr)), 0, 0);

ylims = [floor(min([Ehat_learn2(:); Ehat_orig(:)])/5)*5 0];
xlims = [min(tt_step), max(tt_step)+.025];

lw = 1.5;
subplot(3,1,1) % Head
plot(tt_step, head_curr,'k','Clipping','off','LineWidth',lw); hold on
ylabel('Head')
set(gca,'XColor','w');
ylim(-flip(ylims))
plot(max(xlims)+[0 0], 2+[0 10],'k','LineWidth',1); % 10 deg/s scale bar
% title(sprintf('T_{smooth} %g s', T_smooth));
axis off

subplot(3,1,2) % Eye
plot(tt_step, Ehat_orig, 'Color', learn_color(2,:),'Clipping','off','LineWidth',lw); hold on
plot(tt_step, Ehat_learn2, 'Color', learn_color(3,:),'Clipping','off','LineWidth',lw); hold on
plot(tt_step, Ehat_learn0, 'Color', learn_color(1,:),'Clipping','off','LineWidth',lw); hold on
ylabel('Eye velocity (deg/s')
set(gca,'XColor','w');
ylim(ylims)
plot(T_step+[-0.1 0], -19+[0 0],'k') % Scale bar

axis off

subplot(3,1,3) % PC
plot(tt_step, Phat_orig,  'Color', learn_color(2,:),'Clipping','off','LineWidth',lw); hold on
plot(tt_step, Phat_learn2, 'Color', learn_color(3,:),'Clipping','off','LineWidth',lw); hold on
plot(tt_step, Phat_learn0, 'Color', learn_color(1,:),'Clipping','off','LineWidth',lw); hold on
ylabel('PC (sp/s)')
ylim(ylims)

axis off
linkaxes(findobj(gcf ,'Type','Axes'), 'x')
xlim(xlims)

for jj = 1:3
    subplot(3,1,jj); hold on;
    h = plot(tt_step, tt_step*0,'k--','LineWidth',.5); uistack(h, 'bottom'); 
end

shrink([.4 .7])
fixticks
set(gca,'Fontsize',8)
