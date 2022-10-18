function [hf, hp]  = plotFilters(K)
% function [hf, hp]  = plotFilters(K)

plot_all = 1;
color = 'k';
nCols = 3;

dt = K.dt;

%% 5. Plot resulting filter
hp = [];
hf = figure;

tf = (0:length(K.PH)-1)*dt;
tf_EP = (0:length(K.EP)-1)*dt;
tf_vis = (0:length(K.PR_vel)-1)*dt;
tf_tar_sine = (0:length(K.PT_vel)-1); % TODO: Change based on "enum" option or not

% k_EH
hs = subplot(2, nCols, nCols+1);
hp.EH(1) = plot(tf, K.EH,'k');     hold on

hold on;
h = plot([tf(1) tf(end)], [0 0],'-','Color',.5*[1 1 1]); %uistack(h,'bottom')
ylabel('weight')
title('K_{EH}')
xlabel('Time (s)')

% k_EP
if plot_all
    hs(2) = subplot(2, nCols, nCols+2);
    plot(tf_EP, K.EP,'k');   hold on;
    hp.EP = plot([tf_EP(1) tf_EP(end)], [0 0],'-','Color',.5*[1 1 1]); %uistack(h,'bottom')
    title('K_{EP}')
    xlabel('Time (s)')
end

% k_PH
hs(2,1) = subplot(2, nCols, 1);
hp.PH(1) = plot(tf, K.PH,color);  hold on;
hold on;    h = plot([tf(1) tf(end)], [0 0],'-','Color',.5*[1 1 1]);% uistack(h,'bottom')
ylabel('weight')
H_str = {'H_{vel}'};
legend(hp.PH, H_str(1:length(hp.PH)),'Box','off')
title('K_{PH}');

% k_PR, k_PT
if plot_all
    
    hs(2,2) = subplot(2, nCols, 2); cla
    h = plot([tf_vis(1) tf_vis(end)], [0 0],'-','Color',.5*[1 1 1]); %uistack(h,'bottom')
    hold on;
    hp.PR = plot(tf_vis, K.PR_vel, color);
    hold on; title('Retinal slip (K_{PR})')
    
    % Predictive target sines
    subplot(2, nCols, 3); cla
    hp.PT_sine = [];
    try
    hp.PT_sine(end+1) = plot(tf_tar_sine, K.PT_vel,'-o','Color',0*[.4 0 .6],'MarkerFaceColor',0*[.4 0 .6]); hold on
    hp.PT_sine(end+1) = plot(tf_tar_sine, K.PT_acc*10,'--o','Color',[.4 0 .6]); hold on
    T_sine_str = {'PT_{vel\_sine}','PT_{acc\_sine}'};
    legend(hp.PT_sine, T_sine_str,'Box','off')
    end
    title('Predictive target: sines')
    xlabel('Freq (Hz)')
    
    subplot(2, nCols, 6); cla
    hp.PT_step  = []; T_step_str = {};
    try
        tf_tar_step = (0:length(K.PT_vel_step)-1)*dt;
        hp.PT_step(end+1) = plot(tf_tar_step, K.PT_vel_step,'-','Color',0*[.4 .3 .6]);  hold on
        T_step_str(end+1) = { 'PT_{vel\_step}'};
    end
    try
        tf_tar_step = (0:length(K.PT_acc_step)-1)*dt;
        hp.PT_step(end+1) = plot(tf_tar_step, 100*K.PT_acc_step,'--','Color',[.4 .3 .6]);
        T_step_str(end+1) = {'PT_{acc\_step} * 100'};
    end
    if ~isempty(T_step_str)
        legend(hp.PT_step, T_step_str,'Box','off'); end
    title('Predictive target: steps')
    xlabel('Time (s)')
    
end

drawnow


