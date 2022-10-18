%% Plot PC activity during VORD and cancellation
function [fH_vor_pop, fH_canc_pop, fH_vord, fH_canc, fH_vord_split, fH_canc_split] =...
    plotCancellationPop(I,  K1, K2, K0, pop, B, S)
% [fH_vor_pop, fH_canc_pop, fH_vord, fH_canc, fH_vord_split, fH_canc_split, fH_traces] =...
%     plotCancellationPop(I,  K1, K2, K0, pop, B, S)

%% Plotting parameters
split_target_retina = 1; % plot target and retinal slip separately or not?

color = [.7*[1 1 1];.4*[1 1 1];[0 0 0]]; % Shades of grey
msize = 5;
I.sine_steadystate = 0;

ylims_canc = [-.6 2.8];
% ylims_canc = [0 1.8];
ylims_vord = .8 * [-1 1];

if K1.PF == -1
    ylims_canc = [-3 3];
    ylims_vord = 1.6 * [-1 1];
end

if ~isfield(I,'T_steadystate')
    I.T_steadystate = 6;
end
tt =  0:K1.dt:(I.T_steadystate - K1.dt);

% Function to negate gain if abs(phase) > 90 deg
flipGain = @(gain, ph) gain * ((abs(wrapTo180(ph))<90)*2-1);
learnStr = {'Low', 'Normal','High'};

ncells = pop.ncells;
[P_gain_vor, P_gain_canc,  P_gain_purs] = deal(NaN(3, ncells));
[P_phase_vor, P_phase_canc,  P_phase_purs] = deal(NaN(3, ncells));
[E_gain_vor, E_gain_canc,  E_gain_purs] = deal(NaN(3, 1));
[E_phase_vor, E_phase_canc,  E_phase_purs] = deal(NaN(3, 1));

P_gain_vor_source = [];
P_gain_canc_source = [];
P_phase_canc_source = [];


if pop.scale_total
    
    % Create noise in individual weights
    Es.PH =  (pop.sigmaPH*randn(1, ncells) + 1);
    Es.PR =  (pop.sigmaPR*randn(1,ncells) + 1);
    Es.PE =  (pop.sigmaPE*randn(1,ncells) + 1);
    
    % Make sure mean of all the weights is exactly 1
    Es.PH(end) = ncells - sum(Es.PH(1:end-1));
    Es.PR(end) = ncells - sum(Es.PR(1:end-1));
    Es.PE(end) = ncells - sum(Es.PE(1:end-1));
    
    if ~pop.sigmaPR
        Es.PE = Es.PH;
    end
    
else
    
    % Create noise in individual weights
    Es.PH =  (pop.sigmaPH*randn(nnz(B.bPH_mask), ncells) + 1);
    Es.PR =  (pop.sigmaPR*randn(nnz(B.bPR_vel_mask),ncells) + 1);
    Es.PT_vel = (pop.sigmaPR*randn(nnz(B.bPT_vel_mask),ncells) + 1);
    Es.PT_acc = (pop.sigmaPR*randn(nnz(B.bPT_acc_mask),ncells) + 1);
    Es.PE =  (pop.sigmaPE*randn(nnz(B.bPE_mask),ncells) + 1);
    
    % Make sure mean of all the weights is exactly 1
    Es.PH(:,end) = ncells - sum(Es.PH(:,1:end-1),2);
    Es.PR(:,end) = ncells - sum(Es.PR(:,1:end-1),2);
    Es.PT_vel(:,end) = ncells - sum(Es.PT_vel(:,1:end-1),2);
    Es.PT_acc(:,end) = ncells - sum(Es.PT_acc(:,1:end-1),2);
    Es.PE(:,end) = ncells - sum(Es.PE(:,1:end-1),2);
    
end


Phat_source_canc_mean = [];
Phat_source_canc_mean2 = [];

for jj = 1:3
    
    %% Run for each learning condition
    switch jj
        case 1;         K = K0;
        case 2;         K = K1;
        case 3;         K = K2;
    end
    
    %% Generate population of individual cells
    K_pop = K;
    
    if pop.scale_total
        
        K_pop.PH = K.PH * Es.PH;
        K_pop.PE = K.PE * Es.PE;
        K_pop.PR_vel = K.PR_vel * Es.PR;
        K_pop.PT_vel = K.PT_vel * Es.PR;
        K_pop.PT_acc = K.PT_acc * Es.PR;
        
    else
        
        K_pop.PH = S.scale_H_vel * B.B_PH * bsxfun(@times, K.Kb_eye_pc(B.bPH_mask), Es.PH);
        K_pop.PE = B.B_PE * bsxfun(@times, K.Kb_eye_pc(B.bPE_mask), Es.PE); % Don't use scale here
        K_pop.PR_vel = S.scale_R_vel * B.B_vis * bsxfun(@times, K.Kb_eye_pc(B.bPR_vel_mask), Es.PR);
        K_pop.PT_vel = S.scale_T_vel * bsxfun(@times, K.Kb_eye_pc(B.bPT_vel_mask), Es.PT_vel);
        K_pop.PT_acc = S.scale_T_acc * bsxfun(@times, K.Kb_eye_pc(B.bPT_acc_mask), Es.PT_acc);
        
    end
    
    %% Run model - VORD
    freq = 0.5;

    % Only fit sine waves to the steady state portion of the response
    mask_ss = (length(tt) - round(1/(K.dt*freq)) + 1):length(tt);
        
    head_curr = sin(2*pi*tt*freq);
    target_curr = zeros(size(tt));
    light_on = 0;
    sine_on = freq; % (Hz);
    [Ehat, Phat, Phat_source_vord] = modelClosedloop(K_pop, I, head_curr, target_curr,...
        light_on,  sine_on);
    
    % Fit a sine wave to each PC and eye response
    for cc = 1:ncells
        [gain, P_phase_vor(jj, cc)]  = fitsine(K1.dt, Phat(mask_ss,cc), freq);
        P_gain_vor(jj, cc) = flipGain(gain, P_phase_vor(jj, cc));
    end
    [gain, E_phase_vor(jj)]  = fitsine(K1.dt, Ehat(mask_ss), freq);
    E_gain_vor(jj) = flipGain(gain, E_phase_vor(jj));
    
    % Fit a sine to the average response across all cells
    [P_gain_vor_mean(jj), P_phase_vor_mean(jj)]  = fitsine(K1.dt, ...
        nanmean(Phat(mask_ss,:),2), freq); % Fixed 1/27/19
    if abs(P_phase_vor_mean(jj))>90; P_gain_vor_mean(jj) = - P_gain_vor_mean(jj); end
        
    % Split up by source 
    for kk = 1:size(Phat_source_vord,2)
        Phat_mean = nanmean(Phat_source_vord, 4);
        [gain, ph]  = fitsine(K1.dt, Phat_mean(mask_ss,kk,1), freq); 
        P_gain_vor_source(kk,jj) = flipGain(gain, ph);
    end
    
    %% Run model - cancellation
    head_curr = sin(2*pi*tt*freq);
    target_curr = head_curr;
    light_on = 1;
    [Ehat_canc(:, jj), Phat_canc(:,jj,:), Phat_source_canc] = modelClosedloop(K_pop, I, head_curr, target_curr,...
        light_on,  sine_on);
    
    % Fit a sine wave to each PC and eye response
    for cc = 1:ncells
        [gain, P_phase_canc(jj,cc)]  = fitsine(K1.dt, Phat_canc(mask_ss,jj,cc), freq);
        P_gain_canc(jj, cc) = flipGain(gain, P_phase_canc(jj,cc));
    end
    
    [P_gain_canc_mean(jj), P_phase_canc_mean(jj)]  = fitsine(K1.dt, ...
        nanmean(Phat_canc(mask_ss,jj,:),3), freq);
    P_gain_canc_mean(jj) = flipGain(P_gain_canc_mean(jj), P_phase_canc_mean(jj));
    
    
    [gain, E_phase_canc(jj)]  = fitsine(K1.dt, Ehat_canc(mask_ss,jj), freq);
    E_gain_canc(jj) = flipGain(gain, E_phase_canc(jj));
    

    % Average across cells first
    Phat_source_canc_mean(:,:,:,jj) = nanmean(Phat_source_canc,4);
    
    % Combine across retinal slip and target
    Phat_source_canc_mean2(:,:,:,jj) = nansum(nanmean(Phat_source_canc,4), 3);
    
    if split_target_retina
        for ii = 1:size(Phat_source_canc_mean,2) % Loop over head, visual, and efference
            for ll = 1:size(Phat_source_canc_mean, 3) % Loop over retinal slip and target
                [gain, ph]  = fitsine(K1.dt, Phat_source_canc_mean(mask_ss,ii,ll,jj), freq);
                P_gain_canc_source(ii, jj, ll) = flipGain(gain, ph);
                P_phase_canc_source(ii, jj, ll) = ph;
            end
        end
    else
        
        for ii = 1:size(Phat_source_canc_mean2,2) % Loop over head, visual, and efference
            [gain, ph]  = fitsine(K1.dt, Phat_source_canc_mean2(mask_ss,ii,1,jj), freq);
            P_gain_canc_source(ii, jj, 1) = flipGain(gain, ph);
            P_phase_canc_source(ii, jj, 1) = ph;
            
        end
        
    end
    
    %% Run model - pursuit
    K_pop_purs = K_pop;
    if pop.scale_total
        K_pop_purs.PT_vel = K1.PT_vel * Es.PR;
        K_pop_purs.PT_acc = K1.PT_acc * Es.PR;
    else
        K_pop_purs.PT_vel = S.scale_T_vel * bsxfun(@times, K1.Kb_eye_pc(B.bPT_vel_mask),Es.PT_vel);
        K_pop_purs.PT_acc = S.scale_T_acc * bsxfun(@times, K1.Kb_eye_pc(B.bPT_acc_mask),Es.PT_acc);
    end
    
    target_curr = sin(2*pi*tt*freq);
    head_curr = zeros(size(target_curr));
    light_on = 1;
    [Ehat, Phat, Phat_source_purs] = modelClosedloop(K_pop_purs, I, head_curr, target_curr,...
       light_on,  sine_on);
    
    for cc = 1:ncells
        [gain, P_phase_purs(jj,cc)]  = fitsine(K1.dt, Phat(mask_ss,cc), freq);
        P_gain_purs(jj, cc) = flipGain(gain, P_phase_purs(jj,cc));
    end
    
    [gain, E_phase_purs(jj)]  = fitsine(K1.dt, Ehat(mask_ss), freq);
    E_gain_purs(jj) = flipGain(gain, E_phase_purs(jj)); 
    
end


%% Plot cancellation v. pursuit sens. for PC population
fH_canc_pop = figure;

% Get PC pursuit sensitivity by dividing by mean eye velocity during % pursuit
P_sens_canc = P_gain_canc; % (sp/s per deg/s head) Since the gain of the head is 1
P_sens_purs = bsxfun(@rdivide, P_gain_purs, E_gain_purs); % (sp/s per deg/s eye)

xlims = [0 4];
x_slope = [];
h = [];
for jj = 1:3 %  low norm high plot order
    a = P_sens_purs(jj,:);
    b = P_sens_canc(jj,:);
    % Plot lines of best fit
    x_slope(:,jj) = [a(:) ones(size(a(:)))]\b(:); % A*x = B
    hold on;
    plot(xlims,[xlims(:) ones(2,1)]*x_slope(:,jj),'-','Color',color(jj,:))
    
    % Plot cells
    h(jj) = scatter(P_sens_purs(jj,:), P_sens_canc(jj,:),msize^2, color(jj,:));
end

xlabel('"Eye sensitivity"')
ylabel('"Head sensitivity"')
xlim(xlims); ylim(xlims);
set(gca,'XTick',0:1:xlims(2), 'YTick',0:1:xlims(2))
axis square;
shrink(.4);
fixticks;

%% 1a. Plot VORD modulation of PC firing - individual cells
fH_vor_pop = figure; hold on;
color = I.learn_color;

% Plot dashed line
line([0.5 3.5], [0 0],'Color','k','LineStyle','--'); hold on;

% Loop over VOR gains - low-med-high
for jj = 1:3
    Y = P_gain_vor(jj,:);
    
    % Individual cells
    h = scatter(jj*ones(ncells,1)+randn(ncells,1)*.05, Y, msize^2,...
        'MarkerEdgeColor',color(jj,:), 'MarkerFaceColor',color(jj,:));
    set(h,'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
    ylim([-1 1])
    xlim([.2 3.4]);
    
    ylabel('PC VORD (sp/s per \circ/s)','Interpreter','Tex')
    set(gca,'XTick',1:3,'XTickLabel',learnStr)
end
axis square;
shrink(.4); fixticks; 

%% 2a. Plot VORD modulation of PC firing
fH_vord = figure;

color = I.learn_color;
c_edge_curr = zeros(3,3);
for jj = 1:3

        bar(jj, P_gain_vor_mean(jj), .8,'EdgeColor','none',...% ,c_edge_curr(jj,:)
        'FaceColor',color(jj,:),'LineWidth',2); hold on;
    
end
errorb(1:3, P_gain_vor_mean, sem(P_gain_vor,2));  % NOTE: error bars are sem

ylim(ylims_vord)
set(gca,'YTick',min(ylims_vord):.2:max(ylims_vord))
xlim([.2 3.8]);
% ylabel('PC VORD (sp/s per deg/s)','Interpreter','Tex')
set(gca,'XTick',1:3,'XTickLabel',learnStr)
axis square
shrink(.4); fixticks; 


%% 2b. Plot cancellation modulation, or slope of relationship

% Option 1: Plot slope of relationship between CANCELLATION and pursuit sensitivity 
% P_plot = x_slope(1,:)';

% Option 2: Plot CANCELLATION modulation alone 
P_plot = P_gain_canc_mean(:);%P_gain_canc;

fH_canc = figure; hold on;
color = ones(3,3);
c_edge_curr = I.learn_color;


for jj = 1:3
    bar(jj, mean(P_plot(jj,:)),.8,'EdgeColor',c_edge_curr(jj,:),...
        'FaceColor',color(jj,:),'LineWidth',2); hold on;
end

ylim(ylims_canc)
xlim([.2 3.8]);
set(gca,'YTick',min(ylims_canc):.5:max(ylims_canc))
ylabel('PC canc. (sp/s per deg/s)','Interpreter','Tex')
set(gca,'XTick',1:3,'XTickLabel',learnStr)
axis square
shrink(.4); fixticks;


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
set(gca,'XColor','w');
% ylabel('Contribution to PC VORD')
set(gca,'YTick',min(ylims_vord):.2:max(ylims_vord));
box off; drawnow
set(gca,'XLim',xlims,'YLim',ylims_vord)
shrink([.8 .4]); fixticks;
% offsetAxes;
drawnow;

%%  3b. Plot contributions of Head, Ret slip, and Eye feedback to PC modulation during cancellation - COMBINED TARGET AND RS
%{
fH_canc_split = figure;

color = ones(3,3,3);
c_edge_curr = c_source;
labels_gain = {'L','N','H'};
labels_type = {'Head','Visual','Efference'};
for ii = 1:3 % Loop of head, visual, efference
    for jj = 1:3 % Loop of Low Normal High gain
        
        % Combined target and RS
        yy = sum(P_gain_canc_source(ii, jj, :)); % Not stacked
        x = xs(ii,jj) * ones(size(yy));
        h = bar(x, yy,.8,...
            'EdgeColor',c_edge_curr(jj,:,ii),'FaceColor',color(jj,:,ii),...
            'Clipping','off','LineWidth',2); hold on
        
        if length(h)>1
            set(h(2),'FaceColor',get(h(1),'FaceColor')*.5);
        end
    end
    
    % Plot labels
    text(x-1, ylims_canc(2)*.8,labels_type{ii},'Hor','center');
    
end

% plot(xlims, [0 0], '-k','LineWidth',.5)
set(gca,'YTick',min(ylims_canc):.2:max(ylims_canc))
set(gca,'XColor','w');
ylabel('Contribution to PC canc.')
box off; drawnow
set(gca,'XLim',xlims,'YLim',ylims_canc)

shrink([.8 .4]); fixticks; 
%}


%%  3b. Plot contributions of Head, Ret slip, and Eye feedback to PC modulation during cancellation - PLOT TARGET AND RS CONTRIBUTIONS

fH_canc_split = figure;

color = ones(3,3,3);
c_edge_curr = c_source;
labels_gain = {'L','N','H'};
labels_type = {'Head','Visual','Efference'};
for ii = 1:3 % Loop of head, visual, efference
    for jj = 1:3 % Loop of Low Normal High gain
        
        % Combined target and RS
        yy = squeeze(P_gain_canc_source(ii, jj, :))'; % STACKED
        yy(yy==0) = [];
        h = bar([xs(ii,jj); NaN], [yy; NaN(size(yy))],.8,'Stacked',...
            'EdgeColor',c_edge_curr(jj,:,ii),'FaceColor',color(jj,:,ii),...
            'Clipping','off','LineWidth',2); hold on
        
        if length(h)>1
            set(h(1),'FaceColor',get(h(1),'FaceColor')*.5);
        end
    end
    
    % Plot labels
    text(mean(xs(ii,:)), ylims_canc(2)*.8,labels_type{ii},'Hor','center');
    
end

% plot(xlims, [0 0], '-k','LineWidth',.5)
set(gca,'YTick',min(ylims_canc):.2:max(ylims_canc))
set(gca,'XColor','w');
ylabel('Contribution to PC canc.')
box off; drawnow
set(gca,'XLim',xlims,'YLim',ylims_canc)

shrink([.8 .4]); fixticks; 

%%
%% Plot the traces for each input separately
fH_traces = figure;
colors = [.8 .8 .8; .5 .5 .5; .2 .2 .2];
xt = -1;
if K.PF == 1
    scale_small = [50 10 5];
else
    scale_small = [10  0.5 10];
end
for jj = 1:3  % Loop of Low Normal High gain
    
    text(xt,0, 'Head','Color','r')
    plot(tt, Phat_source_canc_mean(:,1,1,jj), 'Color', colors(jj,:) ); hold on;
    
    text(xt,-2.3, sprintf('Retinal slip x %g', scale_small(1)),'Color','r')
    plot(tt, -2+scale_small(1)*Phat_source_canc_mean(:,2,1,jj),'Color',colors(jj,:)); hold on;
    
    text(xt,-4, sprintf('Target x %g', scale_small(2)),'Color','r')
    plot(tt, -4+scale_small(2)*Phat_source_canc_mean(:,2,2,jj),'Color',colors(jj,:)); hold on;
    
    text(xt,-6.32, sprintf('Efference copy x %g',scale_small(3)),'Color','r')
    plot(tt, -6+scale_small(3)*Phat_source_canc_mean(:,3,1,jj),'Color',colors(jj,:)); hold on;
    
    
    text(xt,-10, 'Purkinje total','Color','r')
    Phat_total = Phat_source_canc_mean(:,1,1,jj) + Phat_source_canc_mean(:,2,1,jj) + Phat_source_canc_mean(:,2,2,jj) + Phat_source_canc_mean(:,3,1,jj);
    plot(tt, -10+Phat_total,'Color',colors(jj,:),'LineWidth',2); hold on;
%     plot(tt,-10+nanmean(Phat_canc(:,jj,:),3),'--','Color',colors(jj,:),'LineWidth',4); % These should be the same
    
    text(xt,-14, 'Eye total x 10','Color','r')
    plot(tt, -14+Ehat_canc(:,jj)*10,'Color',colors(jj,:),'LineWidth',2); hold on;
    
end

title(sprintf('Purkinje cell during cancellation: Positive feedback = %g', K.PF))
fixticks



%}

% %% Plot cancellation v. pursuit sens. for PC population - CHECK HGVPs
% % Only include cells that are HGVPs - less than 45 deg pursuit and canc
% % response, at least 0.33 sens (relative to stimulus) during each
% fH_canc_pop_mask = figure;
% 
% % Get PC pursuit sensitivity by dividing by mean eye velocity during
% % pursuit
% P_sens_canc = P_gain_canc; % (sp/s per deg/s head) Since the gain of the head is 1
% P_sens_purs = bsxfun(@rdivide, P_gain_purs, E_gain_purs); % (sp/s per deg/s eye)
% 
% % Mask based on "Normal" condition only
% hgvp_mask = abs(P_phase_canc(2,:))<=45 & abs(P_phase_purs(2,:))<=45 & P_gain_canc(2,:) > 0.33 & P_gain_purs(2,:) > 0.33;
% 
% xlims = [0 4];
% x_slope = [];
% h = [];
% for jj = 1:3 %  low norm high plot order
%     a = P_sens_purs(jj,hgvp_mask);
%     b = P_sens_canc(jj,hgvp_mask);
%     % Plot lines of best fit
%     x_slope(:,jj) = [a(:) ones(size(a(:)))]\b(:); % A*x = B
%     hold on;
%     plot(xlims,[xlims(:) ones(2,1)]*x_slope(:,jj),'-','Color',color(jj,:))
%     
%     % Plot cells
%     h(jj) = scatter(P_sens_purs(jj,hgvp_mask), P_sens_canc(jj,hgvp_mask),msize^2, color(jj,:));
% end
% 
% xlabel('"Eye sensitivity"')
% ylabel('"Head sensitivity"')
% xlim(xlims); ylim(xlims);
% set(gca,'XTick',0:1:xlims(2), 'YTick',0:1:xlims(2))
% axis square;
% shrink(.4);
% fixticks;

% return

%% DEBUG: plot scale factor for kPR v. pursuit response
%{
figure;

if pop.sigmaPR~=0
    noise_plot = nanmean(Es.PR,1);
else
    noise_plot = Es.PE;
end
% Plot gains
subplot(221)
for jj = 1:3
    plot(noise_plot, P_gain_canc(jj,:),'o','Color',color(jj,:)); hold on
end
xlabel('Noise'); ylabel('Canc gain'); axis square

subplot(222)
for jj = 1:3
    plot(noise_plot, P_gain_purs(jj,:),'o','Color',color(jj,:)); hold on
end
xlabel('Noise'); ylabel('Pursuit gain'); axis square

% Plot phases of PC response
subplot(223)
for jj = 1:3
    plot(noise_plot, P_phase_canc(jj,:),'o','Color',color(jj,:)); hold on
end
plot([min(noise_plot) max(noise_plot)], 45*[1 1],'--k')
plot([min(noise_plot) max(noise_plot)],-45*[1 1],'--k')
xlabel('Noise'); ylabel('Canc phase'); axis square

subplot(224)
for jj = 1:3
    plot(noise_plot, P_phase_purs(jj,:),'o','Color',color(jj,:)); hold on
end

plot([min(noise_plot) max(noise_plot)], 45*[1 1],'--k')
plot([min(noise_plot) max(noise_plot)],-45*[1 1],'--k')
xlabel('Noise'); ylabel('Pursuit phase'); axis square
xlim([min(noise_plot) max(noise_plot)])
linkaxes(get(gcf,'Children'),'x')
fixticks

%% Plot gain v. phase of individual cells
figure;
for jj = 2%1:3
    % plot(P_phase_canc(jj,:), P_gain_canc(jj,:),'ok')
    subplot(121);
    plot(P_phase_purs(jj,:), P_gain_purs(jj,:),'o','Color',color(jj,:)); hold on;
    xlabel('PC phase during pursuit (deg)')
    ylabel('PC gain during pursuit');
    set(gca,'YTick',0:1:3)
    set(gca,'XTick',-90:30:90)
    
    
    subplot(122)
    plot(P_phase_canc(jj,:), P_gain_canc(jj,:),'o','Color',color(jj,:)); hold on;
    
    xlabel('PC phase during cancellation (deg)')
    ylabel('PC gain during cancellation');
    set(gca,'YTick',0:1:3)
    set(gca,'XTick',-90:30:90)
    
end

shrink([.7 .5])
fixticks
linkaxes
xlim([-90 90])
ylim([0 3])
%}

%% Plot VORD phase - EYE
%{
fH_vor_polar = figure;
ax = polaraxes;
for jj=3:-1:1
    %     h = polarplot(deg2rad(P_phase_vor(jj,:)), abs(P_gain_vor(jj,:)),'o','MarkerSize',5); % PC
    
    temp = max(abs(E_gain_vor))*sind(45);
    % h = compass(temp,temp); set(h,'Color','w','LineWidth',4); hold on
    %     h = compass(abs(E_gain_vor(jj))*cosd(E_phase_vor(jj)), abs(E_gain_vor(jj))*sind(E_phase_vor(jj))); % EYE
    
    h = polarplot([0 deg2rad(E_phase_vor(jj))], [0 abs(E_gain_vor(jj))],'-'); hold on% EYE
    
    h(end+1) = polarplot([deg2rad(E_phase_vor(jj))], [abs(E_gain_vor(jj))],'o'); % EYE
    
    set(h,'Color',color(jj,:),'MarkerFaceColor','w','LineWidth',2,'MarkerSize',8); hold on
end
shrink
ax.RLim = [0 1.6];
%}





