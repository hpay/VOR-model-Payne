%% Plot PC activity during VORD and cancellation
function [fH_vord, fH_canc, fH_vord_split, fH_canc_split] = ...
    plotCancellation(I,  K1, K2, K0)
ylims_canc = [-.6 2.2];
ylims_vord = .6 * [-1 1];

tt = (1:8000)*K1.dt-K1.dt; % Allow two cycles

% Function to negate gain if abs(phase) > 90 deg
flipGain = @(gain, ph) gain * ((abs(wrapTo180(ph))<90)*2-1);
learnStr = {'Low', 'Normal','High'};


P_gain_vor = [];
E_gain_vor = [];
P_gain_canc = [];
E_gain_canc = [];

P_gain_vor_source = [];
P_gain_canc_source = [];

for jj = 1:3
    
    switch jj
        case 1;         K = K0;
        case 2;         K = K1;
        case 3;         K = K2;
    end
    
    % Run model - VORD
    head_curr = sin(2*pi*tt*0.5);
    target_curr = zeros(size(tt));
    light_on = 0;
    sine_on = 0.5; % (Hz);
    [Ehat, Phat, Phat_source] = modelClosedloop(K1.dt, head_curr, target_curr,...
        K, I,light_on,  sine_on);
    
    Ehat = Ehat(length(Ehat)/2+1:end);
    Phat = Phat(length(Phat)/2+1:end);
    [P_gain, P_ph]  = fitsine(K1.dt, Phat, 0.5);
    P_gain_vor(jj) = flipGain(P_gain, P_ph);
    [E_gain, E_ph]  = fitsine(K1.dt, Ehat, 0.5);
    E_gain_vor(jj) = flipGain(E_gain, E_ph);
    
    % Split up by source
    for kk = 1:size(Phat_source,2)
        [gain, ph]  = fitsine(K1.dt, Phat_source(round(length(Phat_source)/2):end,kk,1), 0.5);
        P_gain_vor_source(kk,jj) = flipGain(gain, ph);
    end
    
    % Run model - cancellation
    head_curr = sin(2*pi*tt*0.5);
    target_curr = head_curr;
    light_on = 1;
    [Ehat, Phat, Phat_source] = modelClosedloop(K1.dt, head_curr, target_curr,...
        K, I,light_on,  sine_on);
    Ehat = Ehat(length(Ehat)/2+1:end);
    Phat = Phat(length(Phat)/2+1:end);
    [gain, ph] = fitsine(K1.dt, Phat, 0.5);
    P_gain_canc(jj) = flipGain(gain, ph);
    [gain, ph]  = fitsine(K1.dt, Ehat, 0.5);
    E_gain_canc(jj) = flipGain(gain, ph);
    
    % Split up by source
    for kk = 1:size(Phat_source,2)
        for ll = 1:size(Phat_source, 3)
            [gain, ph]  = fitsine(K1.dt, Phat_source(round(length(Phat_source)/2):end,kk,ll), 0.5);
            P_gain_canc_source(kk,jj, ll) = flipGain(gain, ph);
        end
    end
    
end

%% 2a. Plot VORD modulation of PC firing
fH_vord = figure;
I.learn_color = [.7*[1 1 1];.4*[1 1 1];[0 0 0]];

c_curr = I.learn_color;
c_edge_curr = zeros(3,3);
for jj = 1:3
    bar(jj, P_gain_vor(jj), .8,'EdgeColor','none',...% ,c_edge_curr(jj,:)
        'FaceColor',c_curr(jj,:),'LineWidth',2); hold on;
end
ylim(ylims_vord)
set(gca,'YTick',min(ylims_vord):.2:max(ylims_vord))
xlim([.2 3.8]);
% ylabel('PC VORD (sp/s per deg/s)','Interpreter','Tex')
set(gca,'XTick',1:3,'XTickLabel',learnStr)
axis square
shrink(.4); fixticks; cleanticks('y')


%% 2b. Plot CANCELLATION modulation of PC firing
fH_canc = figure; hold on;
c_curr = ones(3,3);
c_edge_curr = I.learn_color;

for jj = 1:3
    bar(jj, P_gain_canc(jj),.8,'EdgeColor',c_edge_curr(jj,:),...
        'FaceColor',c_curr(jj,:),'LineWidth',2); hold on;
end

ylim(ylims_canc)
xlim([.2 3.8]);
set(gca,'YTick',min(ylims_canc):.2:max(ylims_canc))
% ylabel('PC canc. (sp/s per deg/s)','Interpreter','Tex')
set(gca,'XTick',1:3,'XTickLabel',learnStr)
axis square
shrink(.4); fixticks; cleanticks('y')


%%  3a. Plot contributions of Head, Ret slip, and Eye feedback to PC modulation during VORD
fH_vord_split = figure;

c_source = cat(3, I.c(1,:), [190 30 45]/255, I.c(end,:));
c_source =[lighten(c_source, .7); lighten(c_source, .4); c_source];
c_edge_curr = zeros(3,3,3);

dx = 3.5;
xs = (1:3)'*dx*[1 1 1] + [1 1 1]'*(1:3);

for ii = 1:3 % which input (row)
for jj = 1:3 % which learning cond (col)
   
    bar(xs(ii,jj), P_gain_vor_source(ii,jj,1), .8, ...
        'EdgeColor','none','FaceColor',c_source(jj,:,ii),... % c_edge_curr(jj,:,ii)
        'Clipping','off','LineWidth',2); hold on;
end
end
xlims = [3.5 14.5];
% plot(xlims, [0 0], '-k','LineWidth',.5)
set(gca,'XColor','w');
% ylabel('Contribution to PC VORD')
set(gca,'YTick',min(ylims_vord):.2:max(ylims_vord));
box off; drawnow
set(gca,'XLim',xlims,'YLim',ylims_vord)
shrink([.8 .4]); fixticks;
offsetAxes; drawnow; cleanticks('y')

%%  3b. Plot contributions of Head, Ret slip, and Eye feedback to PC modulation during cancellation
fH_canc_split = figure;

c_curr = ones(3,3,3);
c_edge_curr = c_source;


ylims = [-2 2];
for ii = 1:3
    for jj = 1:3 % Loop of Low Normal High gain
%         yy = squeeze(P_gain_canc_source(ii, jj, :)); % Stacked
        yy = sum(P_gain_canc_source(ii, jj, :)); % Not stacked
        x = xs(ii,jj) * ones(size(yy));     

        % Stacked bar for retinal slip (first plane of 3rd dim) and target (second plane)
%         h = bar(x, yy,.8,'stacked',...
%         'EdgeColor',c_edge_curr(jj,:,ii),'FaceColor',c_source(jj,:,ii),...
%         'Clipping','off'); hold on
        
% Normal bar (include target and ret slip)
     h = bar(x, yy,.8,...
        'EdgeColor',c_edge_curr(jj,:,ii),'FaceColor',c_curr(jj,:,ii),...
        'Clipping','off','LineWidth',2); hold on
        
        if length(h)>12
            set(h(2),'FaceColor',get(h(1),'FaceColor')*.5);
        end
    end
end

% plot(xlims, [0 0], '-k','LineWidth',.5)
set(gca,'YTick',min(ylims_canc):.2:max(ylims_canc))
set(gca,'XColor','w'); 
% ylabel('Contribution to PC canc.')
box off; drawnow
set(gca,'XLim',xlims,'YLim',ylims_canc)

shrink([.8 .4]); fixticks; cleanticks('y')


%% TODO:
% Plot predictions for model in single trial differences in cancellation performance




