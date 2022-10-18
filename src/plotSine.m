function plotSine(K1, K2, K0, I, learn_color)
% (temmporarly?) helper function -- plot VORD sine for PC

freq = 0.5;
dt = K1.dt;
tt_sine = (0:2/dt)'*dt;
head_curr = sin(2*pi*freq*tt_sine);

% TODO: make sure these match!!!
[Ehat_orig, Phat_orig] = modelImpulse(getImpulse(K1, I, 1), I, head_curr, zeros(size(head_curr)), 0, freq);
[Ehat_learn2, Phat_learn2] = modelImpulse(getImpulse(K2, I, 1), I, head_curr, zeros(size(head_curr)), 0, freq);
[Ehat_learn0, Phat_learn0] = modelImpulse(getImpulse(K0, I, 1), I, head_curr, zeros(size(head_curr)), 0, freq);


% [Ehat_orig, Phat_orig] = modelClosedloop(K1.dt, head_curr, zeros(size(head_curr)), ...
%     K1, I, 0, freq);
% [Ehat_learn2, Phat_learn2] = modelClosedloop(K1.dt, head_curr, zeros(size(head_curr)), ...
%     K2, I, 0, freq);
% [Ehat_learn0, Phat_learn0] = modelClosedloop(K1.dt, head_curr, zeros(size(head_curr)), ...
%     K0, I, 0, freq);
% tt_sine2 = (0:4/dt+1)'*dt;

ylims = [-1.5 1.5];
subplot(3,1,1) % Head
plot(tt_sine, head_curr,'k','Clipping','off')
ylabel('Head')
set(gca,'XColor','w');
ylim(ylims)

subplot(3,1,2) % Eye
plot(tt_sine, Ehat_orig, 'Color', learn_color(2,:),'Clipping','off'); hold on
plot(tt_sine, Ehat_learn2, 'Color', learn_color(3,:),'Clipping','off'); hold on
plot(tt_sine, Ehat_learn0, 'Color', learn_color(1,:),'Clipping','off'); hold on
ylabel('Eye velocity (deg/s')
set(gca,'XColor','w');
ylim(ylims)

subplot(3,1,3) % PC
plot(tt_sine, Phat_orig,  'Color', learn_color(2,:),'Clipping','off'); hold on
plot(tt_sine, Phat_learn2, 'Color', learn_color(3,:),'Clipping','off'); hold on
plot(tt_sine, Phat_learn0, 'Color', learn_color(1,:),'Clipping','off'); hold on
ylabel('PC (sp/s)')
ylim(ylims)

linkaxes(findobj(gcf ,'Type','Axes'), 'x')
xlim([tt_sine(1) tt_sine(end)]);


for jj = 1:3
    subplot(3,1,jj); hold on; h = plot(tt_sine, tt_sine*0,'--k'); uistack(h, 'bottom'); end
shrink([.4 .7])
fixticks
set(gca,'Fontsize',8)