function PLOT_MODEL(filepathname, option_final)
% Plot summary results from models fit (using RUN_MODEL.m) using multiple
% levels of positive feedback
%
% Provide full path to folder with results
%
% Hannah Payne, 2019
% Raymond and Goldman Labs
close all
ptick = 0.5; % scaling for ticks, increase to make longer
if ~exist('option_final','var') 
    option_final = 0;
end
if option_final
    paperDefaults
end
set(groot,'defaultAxesXTickLabelRotationMode','manual')
set(groot,'defaultAxesYTickLabelRotationMode','manual')
set(groot,'defaultAxesZTickLabelRotationMode','manual')

% Specific model result to analyze
[~, filename] = fileparts(filepathname);

% Add helper function folders to search path
code_folder = fileparts(fileparts(which(mfilename)));
addpath(genpath(fullfile(code_folder, 'src')))
folder_root = fullfile(code_folder,'results');

% Load the full model
D = load(fullfile(folder_root,filename,[filename,'.mat']));
I = D.I; % Model info/parameters
B = D.B; % Model basis set
E = D.E; % Model fit errors
S = D.S; % Scale factors for normalizing data going into model fit
R = D.R; % Regularization

% Change colors
I.c = flipud(turbo(12));
I.c(end,:) = [];

% May need to fix path in parameters file (e.g. if current pathname doesn't match original pathname)
I.figures_path = fullfile(folder_root, filename);
if option_final
    dropbox_root = I.figures_path(1:strfind(I.figures_path, 'Dropbox')+6);
    I.figures_path = fullfile(dropbox_root,'rlab','model','Figures');
end
I.impulse_or_closedloop = 0; 

% data_path = fullfile(folder_root, 'VOR-model','data');
% [data_path_temp1, data_path_temp2] = fileparts(data_path);
% data_path = fullfile(fileparts(which(mfilename)), 'data');
data_path = fullfile(code_folder, 'data');

[conds, nConds_JR, tts, head, target, hevel, PC, sines, light, ...
    dt, n_cells, RR_data] = loadJR_RR_combined(I, data_path);
[conds, tts, head, target, hevel, PC, sines, light] = addStep(I, conds, ...
    tts, head, target, hevel, PC, sines, light);

% Select model fits to plot: fine-tuned visual weights
K1 = D.K1;
K0 = D.K0;
K2 = D.K2;

% Positive feedback weights
PFs = I.PFs;

%% Plot weights before learning
pos_weights = [10 10 8 8]; %pos_weights = [10 10 8 9];
figure('Units','centi','Pos',pos_weights); 

% Set xlims
maxt = I.max_time_basis*1000;      % (ms)
tt_max_R = I.max_time_basis_vis*1000;     % (ms)
tt_max_E = I.max_time_basis_EP*1000;     % (ms)
tt_max_step = 200;  % (ms)

% Set ylims
ymax_PH = 0.25;
ymax_EH = 0.8;
ymax_PR = 0.006;
ymax_EP = 0.04;
ylim_PRstep = [-0.2 0.5];
ylim_PHstep = [-.7 1.5];
ylim_Estep = [-.7 1.5];
dPH = 0.125;
dEH = 0.4;
dPR = 0.003;
dEP = .02;


% Create time vectors
ii = 1;
tt_H = (1:length(K1(ii).PH))*dt*1000;
tt_R = (1:length(K1(ii).PR_vel))*dt*1000;
tt_E = (1:length(K1(ii).EP))*dt*1000;

% Timing for step
tbefore = 25;                                           % (ms)
tt_step = (0:1000*dt:(tt_max_step+tbefore))  - tbefore; % (ms)
step_offset = .7; % For plotting example step input
step_scale = 1;

% Create step input
step_in = zeros(1,length(tt_step)); step_in(round(tbefore/1000/dt):end) = 1;
step_in = smooth(step_in,round(.01/dt)); % 10 ms smooth
step_in = smooth(step_in,round(.01/dt)); % 10 ms smooth
hs = gobjects(4,2);
lw=.5;
for ii = flip(find(ismember(PFs(:),(0:.1:1))))'
    if length(K1)<ii || isempty(K1(ii).PH); continue; end

    
  
    si = 0;
     % kEH linear filter
    
    si = si+1;
    hs(si,1) = subplot(4,2,si*2-1);
    plot(tt_H(tt_H<maxt), -K1(ii).EH(tt_H<maxt),'LineWidth',lw,'Color',I.c(ii,:),'Clipping','off'); hold on
    set(gca,'XLim',[tt_H(1) maxt])    
%         text(18, .65, 'No feedback','Color', I.c(1,:), 'FontSize',7)
%     text(18, .3, 'Strong feedback','Color', I.c(end,:), 'FontSize',7)
         
    % kEH step response
    hs(si,2) = subplot(4,2,si*2);
    hold on;
    step_curr = filter(-K1(ii).EH, 1, step_in);          % Step response
    plot(tt_step, step_curr,'Color',I.c(ii,:),'Clipping','off','LineWidth',lw); hold on;
    plot(tt_step, step_offset+step_scale*step_in,'k','Clipping','off');     hold on; % Step input

    
    % kEP linear filter
    si = si+1;
    hs(si,1) = subplot(4,2,si*2-1);
    plot(tt_E(tt_E<tt_max_E), K1(ii).EP(tt_E<tt_max_E),'LineWidth',lw,'Color',I.c(ii,:),'Clipping','off'); hold on
    set(gca,'XLim',[tt_E(1) tt_max_E])

    % kEP step response
    hs(si,2) = subplot(4,2,si*2);
    hold on;
    step_curr = filter(K1(ii).EP, 1, step_in);          % Step response
    plot(tt_step, step_curr,'LineWidth',lw,'Color',I.c(ii,:),'Clipping','off');
   

    
    % kPH linear filter
    si = si+1;
    hs(si,1) = subplot(4,2,si*2-1);
    plot(tt_H(tt_H<maxt), K1(ii).PH(tt_H<maxt),'LineWidth',lw,'Color',I.c(ii,:),'Clipping','off'); hold on
    set(gca,'XLim',[tt_H(1) maxt])
    
    % kPH step response
    hs(si,2) = subplot(4,2,si*2); hold on;
    step_curr = filter(K1(ii).PH, 1, step_in);          % Step response
    plot(tt_step, step_curr,'LineWidth',lw,'Color',I.c(ii,:),'Clipping','off');
    
    % kPR linear filter
    si = si+1;
    hs(si,1) = subplot(4,2,si*2-1);
    plot(tt_R(tt_R<tt_max_R), K1(ii).PR_vel(tt_R<tt_max_R),'LineWidth',lw, 'Color',I.c(ii,:),'Clipping','off'); hold on
    set(gca,'XLim',[tt_R(1) tt_max_R])
    
    % kPR step response
    hs(si,2) = subplot(4,2,si*2);
    hold on;
    step_curr = filter(K1(ii).PR_vel, 1, step_in);          % Step response
    plot(tt_step, step_curr,'LineWidth',lw,'Color',I.c(ii,:),'Clipping','off');
    
end

for ii = 1:numel(hs)
    hs(ii).XRuler.TickLabelGapOffset = -1;
    hs(ii).YRuler.TickLabelGapOffset = 0;
    if mod(ii,4)==0
        xlabel(hs(ii), 'Time (ms)');
    end
    if ii<5
        text(hs(ii), -.34, .5 , 'Filter strength','Units','norm','Rot',90,'Horiz','center', 'FontSize',7)
    end
end

% set(hs(1:2,1), 'XLim',[0 maxt],''XTick',0:10:maxt);
% set(hs(3,1), 'XLim',[0 tt_max_R],'XTick',0:100:tt_max_R);
% set(hs(4,1), 'XLim',[0 tt_max_E],'XTick',0:50:tt_max_E);
set(hs(:,2), 'XLim',[-tbefore tt_max_step], 'XTick',0:50:tt_max_step);
set(hs(:),'color','none')
set(hs(3,1), 'YLim',[-ymax_PH ymax_PH],'YTick', -ymax_PH:dPH:ymax_PH);
set(hs(1,1), 'YLim',[-ymax_EH ymax_EH],'YTick', -ymax_EH:dEH:ymax_EH);
set(hs(4,1), 'YLim',[-ymax_PR ymax_PR], 'YTick', -ymax_PR:dPR:ymax_PR);
set(hs(2,1), 'YLim',[-ymax_EP ymax_EP], 'YTick', -ymax_EP:dEP:ymax_EP);
set(hs(1:2,2), 'YLim',ylim_Estep);
set(hs(3,2), 'YLim',ylim_PHstep);
set(hs(4,2),'YLim',ylim_PRstep)
% set(hs(2,2),'YLim',[-.3 ymax_Estep])

if option_final
    fixticks(.5)
    cleanticks
    export_fig(fullfile(I.figures_path,'baseW.pdf'),'-nocrop','-painters')
else
    ylabel(hs(1,1), 'k_{EH}');
    ylabel(hs(2,1), 'k_{EP}');
    ylabel(hs(3,1), 'k_{PH}');
    ylabel(hs(4,1), 'k_{PR}');
    my_export_fig(fullfile(I.figures_path,'baseW.pdf'))
end


%% ============Run prediction from stimulation to PC==============
%% ============Run prediction from stimulation to PC==============

dt = .0005;
T_total = .350;     % (s) total trace duration
T_stim = .025;      % (s) stim duration
T_start = .1;       % (s) stim start
tt_H = (0:dt:T_total)';
[head_curr, target_curr, pc_stim] = deal(zeros(size(tt_H)));
pc_stim(round(T_start/dt):round((T_start+T_stim)/dt)) = 1; % 25 ms pulse
light_on = 0;
sine_on = 0;


t_space = range(tt_H)*1.1; % spacing from one plot to next;
figure('Units','centi','Pos',[15 10 8 2])
ii_plot = find(ismember(round(PFs(:)'*10),round([0 0.3 0.6 0.9 1]*10))); % Feedback strengths to plot
for ii_temp = 1:length(ii_plot)
    ii = ii_plot(ii_temp);
    
    K = K1(ii);
    if isempty(K.PH); continue; end
    [Ehat, Phat, Phat_source] = modelClosedloop(K, I, head_curr, target_curr,...
        light_on,  sine_on, pc_stim);
    Ehat_norm = Ehat/max(Ehat);
    plot(tt_H + t_space*ii_temp, Ehat_norm, 'Color', I.c(ii,:)); hold on
    h = rectangle('Position',[t_space*ii_temp+T_start  0  T_stim 1],...
        'FaceColor',.8*[1 1 1],'EdgeColor','none');
    uistack(h,'bottom')
end
axis off
xlim([t_space t_space*(ii_temp+1)])

if option_final
    export_fig(fullfile(I.figures_path, 'stimSummary.pdf'))
else
    my_export_fig(fullfile(I.figures_path, 'stimSummary.pdf'))
end

%% Replot final baseline fits 
mask = 1:25; % Choose conditions to plot
hf_basefit =  plotBaselineResultsPair(K1([1 11]), I, conds(mask), tts(mask), head(mask), target(mask), hevel(mask), PC(mask), light(mask), sines(mask));
set(hf_basefit,'Position',[1   1   650   775])
if option_final
    fixticks
    export_fig(hf_basefit, fullfile(I.figures_path, sprintf('baseFit_tuneR.pdf')),'-nocrop');
else    
    my_export_fig(hf_basefit, fullfile(I.figures_path, sprintf('baseFit_tuneR.pdf')));
end




%% Show instability - VOR response with 10% change in weights
iis = find(ismember(PFs,[0 0.9 1]));
for ii = iis(:)'
    if ~isempty(K1(ii).PH)
        figure('Units','centi','Pos',[6   9    4 4]);
        plotInstability(I,  K1(ii), I.c(ii,:));  
        ylim([-25 5])
        xlim(xlim+[0 .05])
        axis off
        if ii==1 
            plot([-.1 0]+max(xlim), -20*[1 1],'k') ;
            plot([0 0]+max(xlim), -20+[0 10],'k')
        end
        if option_final
            fixticks;
            cleanticks;
            cleanticks('x');
            export_fig(fullfile(I.figures_path, sprintf('instability_PF%g.pdf',I.PFs(ii))),'-nocrop')
        else
            my_export_fig(fullfile(I.figures_path, sprintf('instability_PF%g.pdf',I.PFs(ii))))
        end
    end
end


%% Report the gain of visual pursuit when 90% of the flocculus is "lesioned"
for ii = [1 11]
% Make copy of model and set kEP to 20% of normal
K1_curr = K1(ii);
K1_lesion = K1_curr;
K1_lesion.EP = K1_curr.EP*0.2;

% Run simulation
amp = 10;
sine_freq = 0.5;
tt = (0:K1_curr.dt:2);
head_curr = zeros(size(tt));
target_curr = amp*sin(2*pi*sine_freq*tt);

[Ehat, Phat] = modelClosedloop(K1_curr, I, head_curr, target_curr, ...
    1, sine_freq);

% Run lesioned model
[Ehat_lesion, Phat_lesion] = modelClosedloop(K1_lesion, I, head_curr, target_curr, ...
    1, sine_freq);

% Fit a sine to each
[amp_ctrl, phase] =  fitsine(tt, Ehat, sine_freq);
[amp_lesion, phase] =  fitsine(tt, Ehat_lesion, sine_freq);

fprintf('Percent of pursuit gain in intact animals, feedback gain = %.1g: %.2f%%\n', K1_curr.PF,(amp_lesion/amp_ctrl)*100)
end

%% Re-plot the frequency response compared to data
linespec = '-o';
fsize = [10.6627   10.2923    3.3  4.4]; % cm

% Plot freq data
hF_freq = figure('Units','centi','Pos',fsize); drawnow
plotFreq(RR_data.freqs, RR_data.gains1, RR_data.phases1+180, linespec, I.learn_color(2,:), I.learn_color(2,:));
plotFreq(RR_data.freqs, RR_data.gains2, RR_data.phases2+180, linespec, I.learn_color(3,:), I.learn_color(3,:));
% plotFreq(RR_data.freqs, RR_data.gains0, RR_data.phases0+180, linespec, I.learn_color(1,:), I.learn_color(1,:));
xlabel('Frequency (Hz)');   fixticks(ptick); drawnow;
ha = findobj(hF_freq,'Type','Axes'); 
pbaspect(ha(1),[1.5 1 1])
pbaspect(ha(2),[1.5 1 1])

if option_final
    export_fig(hF_freq, fullfile(I.figures_path, 'learn_freq_data.pdf'),'-nocrop');
else
    my_export_fig(hF_freq, fullfile(I.figures_path, 'learn_freq_data.pdf'));
end

for ii = find(ismember(PFs(:)',[0 1])) % Only plot for pos fdbk = 0 or 1
    
    hF_freq = figure('Units','centi','Pos',fsize);
    
    % Get frequency response after learning
    [~, ~, E_gain1, E_phase1, P_gain1, P_phase1] = getFreq(K1(ii), I, RR_data.freqs);
    [~, ~, E_gain2, E_phase2, P_gain2, P_phase2] = getFreq(K2(ii), I, RR_data.freqs);
    [~, ~, E_gain0, E_phase0, P_gain0, P_phase0] = getFreq(K0(ii), I, RR_data.freqs);
    
    % Plot eye frequency response after learning
    plotFreq(RR_data.freqs, E_gain1*RR_data.scaleJR2RR, E_phase1, linespec, I.learn_color(2,:), I.learn_color(2,:));
    plotFreq(RR_data.freqs, E_gain2*RR_data.scaleJR2RR, E_phase2, linespec, I.learn_color(3,:), I.learn_color(3,:));
%     plotFreq(RR_data.freqs, E_gain0*RR_data.scaleJR2RR, E_phase0, linespec, I.learn_color(1,:), I.learn_color(1,:));
    xlabel('Frequency (Hz)');  fixticks(ptick); drawnow;
    ha = findobj(hF_freq,'Type','Axes');
    pbaspect(ha(1),[1.5 1 1])
    pbaspect(ha(2),[1.5 1 1])
    if option_final
        export_fig(hF_freq, fullfile(I.figures_path, sprintf('learn_freq%g.pdf', I.PFs(ii))),'-nocrop');
    else
        my_export_fig(hF_freq, fullfile(I.figures_path, sprintf('learn_freq%g.pdf', I.PFs(ii))));
    end
end




%% Run simulation with visual occlusion/target blink
for ii = find(ismember(PFs(:)',[0 1]))
    figure('Units','centi','Pos',[15 10 6 4])
    plotVisualOcclusion(I,  K1(ii), I.c(ii,:));
    if option_final
        export_fig(fullfile(I.figures_path, sprintf('visualOcclusion_PF%g.pdf',I.PFs(ii))),'-nocrop')
    else
        my_export_fig(fullfile(I.figures_path, sprintf('visualOcclusion_PF%g.pdf',I.PFs(ii))))
    end
end

%% Plot PC activity during VORD and cancellation

% Check if no brainstem plasticity file is present
flag_plot_nobrainstem = 0;
    filepathname_nobrainstem = 'D:\hannah\Dropbox\rlab\model\TO_SHARE\20211001_init_nobrainstem';
% filepathname_nobrainstem = fullfile(folder_root,[filename '_nobrainstem']);
[~, filename_nobrainstem] = fileparts(filepathname_nobrainstem);
if exist(filepathname_nobrainstem,'dir')
    flag_plot_nobrainstem = 1;
    D_nobs = load(fullfile(filepathname_nobrainstem,[filename_nobrainstem '.mat']));    
end

flag_plot_noLTP = 0;
if exist([filepathname,'_noLTP'],'dir')
    flag_plot_noLTP = 1;
    [folder_root, filename] = fileparts(filepathname);
    D_noLTP = load(fullfile(folder_root,[filename '_noLTP'],[filename,'_noLTP.mat']));
    D_noLTD = load(fullfile(folder_root,[filename '_noLTD'],[filename,'_noLTD.mat']));
    
end
for ii = find(ismember(PFs(:)',[0 1]))
    %%
    if ismember(PFs(ii), I.runPFs)
        rng(1)
        % Set sigmas for noise in simulated population of PCS
        if PFs(ii)==1
            pop.sigmaPH = 0.1;
            pop.sigmaPR = 20;
            pop.sigmaPE = 0;
        elseif PFs(ii)==-1
            pop.sigmaPH = 0.1;
            pop.sigmaPR = 0.25;
            pop.sigmaPE = 0;
        else
            pop.sigmaPH = .25;
            pop.sigmaPR = 0.4;
            pop.sigmaPE = 0;
        end
        pop.scale_total = 1; % If 0, randomly scale each basis function individually
        pop.ncells = 20;
        
        pop0 = [];
        pop0.sigmaPH = 0;
        pop0.sigmaPR = 0;
        pop0.sigmaPE = 0;
        pop0.scale_total = 1; % If 0, randomly scale each basis function individually
        pop0.ncells = 1;
        
        % Plot Gain Down and Gain Up separately
        [fH_vord, fH_canc, fH_traces, P_gain_vor_mean, P_gain_canc_mean] = ...
            plotVorCancSplit(I,  K1(ii), K2(ii), K0(ii),  pop, B, S); % option_final
        
        if flag_plot_nobrainstem
            [~, ~, ~, P_gain_vor_mean_nobs, P_gain_canc_mean_nobs] = ...
                plotVorCancSplit(I,  D_nobs.K1(ii), D_nobs.K2(ii), D_nobs.K0(ii),  pop0, B, S); % option_final
            
            % Add dashed line to the gain increase figure for VOR cancellation
            figure(fH_canc(2))
            h=get(fH_canc(2),'children');
            axes(h(2))
            hold on
            plot(1+.4*[-1 1], P_gain_canc_mean_nobs(2)*[1 1],':','Color','r','Linewidth',1)
            plot(2+.4*[-1 1], P_gain_canc_mean_nobs(3)*[1 1],':','Color','r','Linewidth',1)
            
        end
        
        
        types = {'LOW','HIGH'};
        for gg = 1:length(fH_canc)
            type = types{gg};
            if option_final
                export_fig(fH_vord(gg), fullfile(I.figures_path, sprintf('learn_vord_PF%g%s.pdf', PFs(ii), type)),'-painters','-nocrop')
                export_fig(fH_canc(gg), fullfile(I.figures_path, sprintf('learn_canc_PF%g%s.pdf', PFs(ii), type)),'-painters','-nocrop')
                export_fig(fH_traces, fullfile(I.figures_path, sprintf('learn_canc_traces_PF%g.pdf', PFs(ii))),'-painters','-nocrop')
            else
                my_export_fig(fH_vord(gg), fullfile(I.figures_path, sprintf('learn_vord_PF%g%s.pdf', PFs(ii), type)))
                my_export_fig(fH_canc(gg), fullfile(I.figures_path, sprintf('learn_canc_PF%g%s.pdf', PFs(ii), type)))
                my_export_fig(fH_traces, fullfile(I.figures_path, sprintf('learn_canc_traces_PF%g.pdf', PFs(ii))))
            end
        end
        
        %% Run the noLTD and noLTP models
        if flag_plot_noLTP
            flag_gain_up = 1;
            [fH_vord, fH_canc] = plotVorCanc(D_noLTD.I,  D_noLTD.K1(ii), D_noLTD.K2(ii), D_noLTD.K0(ii), flag_gain_up);
            if option_final
                export_fig(fH_vord, fullfile(I.figures_path, sprintf('learn_vord_PF%g_noLTD.pdf', PFs(ii))),'-nocrop')
                export_fig(fH_canc, fullfile(I.figures_path, sprintf('learn_canc_PF%g_noLTD.pdf', PFs(ii))),'-nocrop')
            else
                title(findall(fH_vord,'type','axes'), sprintf('VORD, no LTD, PF = %g', PFs(ii)))
                title(findall(fH_canc,'type','axes'), sprintf('CANC, no LTD, PF = %g', PFs(ii)))
            end
            [fH_vord, fH_canc] = plotVorCanc(D_noLTP.I,  D_noLTP.K1(ii), D_noLTP.K2(ii), D_noLTP.K0(ii), flag_gain_up);
            if option_final
                export_fig(fH_vord, fullfile(I.figures_path, sprintf('learn_vord_PF%g_noLTP.pdf', PFs(ii))),'-nocrop')
                export_fig(fH_canc, fullfile(I.figures_path, sprintf('learn_canc_PF%g_noLTP.pdf', PFs(ii))),'-nocrop')
            else
                
                title(findall(fH_vord,'type','axes'), sprintf('VORD, no LTP, PF = %g', PFs(ii)))
                title(findall(fH_canc,'type','axes'), sprintf('CANC, no LTP, PF = %g', PFs(ii)))
            end
        end
    end
end

%% Plot errors
figure('Units','centimeters','Pos', [10   12    8 2.35]); hold on;
cond_mask = 1:I.nConds_JR;
nrmse_pc = NaN(size(I.PFs));
nrmse_eye = NaN(size(I.PFs));
bw = .06;
for ii = 1:length(PFs)
    hs  = subplot(1,2,1);
    % Use the updated error after fine tuning visual params
    if isempty(E(ii).nrmse_pc_vis); continue; end
    nrmse_pc(ii) = nanmean(E(ii).nrmse_pc_vis(cond_mask));
    nrmse_eye(ii) = nanmean(E(ii).nrmse_eye_vis(cond_mask));
    bar(PFs(ii), nrmse_pc(ii),bw, 'FaceColor',I.c(ii,:),'EdgeColor','none'); hold on;
    
    hs(2) = subplot(1,2,2);
    bar(PFs(ii), nrmse_eye(ii),bw, 'FaceColor',I.c(ii,:),'EdgeColor','none'); hold on;
end

ylims = [0 .05];
linkaxes;
set(hs, 'YTick',ylims(1):ylims(2)/5:ylims(2));
ylim(ylims)
xlim([-.1+min(PFs) max(PFs)+.1])

% Calculate total error
err_total = nrmse_pc + nrmse_eye;
temp = nanmax((err_total-nanmean(err_total))/nanmean(err_total));
fprintf('Max difference in error = %.3g%% from mean\n',temp*100);

% Make axis labels
set(hs,'XTick',sort(PFs))
PF_str = arrayfun(@num2str, sort(PFs),'Uni',false);
PF_str(2:10) = {''};
set(hs,'XTickLabel',PF_str)

hs(1).XRuler.TickLabelGapOffset = -1;
hs(1).YRuler.TickLabelGapOffset = 0;
hs(2).XRuler.TickLabelGapOffset = -1;
hs(2).YRuler.TickLabelGapOffset = 0;

if option_final
     fixticks(ptick);
cleanticks
export_fig(fullfile(I.figures_path,'baseError.pdf'))
else
    ylabel(hs(1), {'Model error:','Purkinje cell rate'})
    ylabel(hs(2),'Eye velocity')
    xlabel(hs(1),'Positive feedback')
    xlabel(hs(2),'Positive feedback')
    my_export_fig(fullfile(I.figures_path,'baseError.pdf'))
end



%% Plot errors for PF = 0 and PF = 1 divided by type
nconds_err = 25;
conds_err = conds(1:nconds_err);
errTemp = NaN(nconds_err,2);
errTemp(:,1, 1) = E(1).nrmse_pc_vis;
errTemp(:,2, 1) = E(11).nrmse_pc_vis;

errTemp(:,1, 2) = E(1).nrmse_eye_vis;
errTemp(:,2, 2) = E(11).nrmse_eye_vis;


xticks = (1:length(conds_err));
indDark = contains(conds_err,'dark');
ind2 = contains(conds_err,'x2');
ind0 = contains(conds_err,'x0');
indPurs = contains(conds_err,'pursuit');
condsOrd = [conds_err(indDark); conds_err(ind2); conds_err(ind0); conds_err(indPurs)];

errTemp = [errTemp(indDark,:,:); errTemp(ind2,:,:); errTemp(ind0,:,:); errTemp(indPurs,:,:)];


xtext = strrep(condsOrd, '_',' ');
xtext = strrep(xtext, 'ms', ' ms');
xtext = strrep(xtext, 'Hz', ' Hz');
xtext = strrep(xtext, '05', '0.5');
xtext = strrep(xtext, 'pursuit', 'purs.');

figure; 

hs = subplot(2,1,1); % PC
h = bar(xticks, errTemp(:,:,1),'Clipping','off');
h(1).FaceColor = I.c(1,:);
h(2).FaceColor = I.c(11,:);
set(h,'EdgeColor', 'none');
set(h,'BarWidth', 1);
set(gca,'XTick',[],'XTickLabel','')
ylabel('Purkinje cell error')
ylim([0 .1])
xlim([min(xticks)-.67 max(xticks)+.3])
set(hs(1),'Pos',get(hs(1),'Pos')+[0 .05 0 0])

hs(2) = subplot(2,1,2); % eye
h = bar(xticks, errTemp(:,:,1),'Clipping','off');
h(1).FaceColor = I.c(1,:);
h(2).FaceColor = I.c(11,:);
set(h,'EdgeColor', 'none');
set(h,'BarWidth', 1);
set(gca,'XTick',[],'XTickLabel','')
for i = 1:length(xticks)
    text( xticks(i),-.01, xtext(i),'FontSize',5,'HorizontalAlignment','right','Rotation',60)
end
ylabel('Eye velocity error')
ylim([0 .1])
xlim([min(xticks)-.67 max(xticks)+.3])
set(hs(2),'Pos',get(hs(2),'Pos')+[0 .1 0 0])
if option_final
    shrink
    fixticks(ptick)
    box off
    export_fig(fullfile(I.figures_path,'baseErrorAll.pdf'))
end





%% Plot change in vestibular weights after learning
tt_H = 1000*((1:length(K1(1).PH))*dt-dt);
ymax1 = 0.12; dy = .05;
pos_learn = [10 10 3.9 5.8];
hs = [];
alpha = 1;
maxt = 50; % ms to plot
dx = 10;   % ticks
types = {'LOW','HIGH'};
for jj =1:length(types)
    
    
    % Normal learning
    K_baseline = K1;
    if strcmp(types{jj},'LOW');     K_learn = K0;
    else                            K_learn = K2;
    end
    
    for ii = find(ismember(PFs(:)', [1 0]))
        figure('Units','centi','Pos',pos_learn)
        if length(K_learn)<ii || isempty(K_learn(ii).PH); continue; end
        hs(1) = subplot(2,1,1);
        %         plot(tt_H(tt_H<=maxt+eps), K_learn(ii).PH(tt_H<=maxt+eps) - K_baseline(ii).PH(tt_H<=maxt),'Color',I.c(ii,:),'Clipping','off'); hold on
        area(tt_H(tt_H<=maxt), K_learn(ii).PH(tt_H<=maxt)-K_baseline(ii).PH(tt_H<=maxt),'FaceColor',I.c(ii,:),'EdgeColor','none','FaceAlpha',alpha,'Clip','on','ShowBaseLine','off'); hold on
        
        hs(2) = subplot(2,1,2);
        %         plot(tt_H(tt_H<=tt_max_H+eps), -(K_learn(ii).EH(tt_H<=tt_max_H+eps) - K_baseline(ii).EH(tt_H<=tt_max_H)),'Color',I.c(ii,:),'Clipping','off'); hold on
        area(tt_H(tt_H<=maxt), -(K_learn(ii).EH(tt_H<=maxt)-K_baseline(ii).EH(tt_H<=maxt)),'FaceColor',I.c(ii,:),'EdgeColor','none','FaceAlpha',alpha,'Clip','on','ShowBaseLine','off'); hold on
        
        axis(hs,'square')
        
        set(hs,'XLim',[0 maxt],'XTick',0:dx:maxt, 'Box','off')
        set(hs(1:2), 'YLim',[-ymax1 ymax1],'YTick',-floor(ymax1/dy)*dy:dy:ymax1)
        subplot(hs(1))
        h = plot(tt_H,zeros(size(tt_H)),'-','Color',0*[1 1 1],'LineWidth',.5);
%         uistack(h,'bottom');
        subplot(hs(2))
        h = plot(tt_H,zeros(size(tt_H)),'-','Color',0*[1 1 1],'LineWidth',.5); 
%         uistack(h,'bottom');
        xlabel('Time (ms)')        
             
        if option_final
            fixticks(ptick)
            cleanticks
            export_fig(fullfile(I.figures_path,sprintf('learnDeltaW_%s_PF%g.pdf', types{jj}, PFs(ii)) ),'-nocrop')
        else
            title(hs(1), types{jj})
            ylabel(hs(1), {'\DeltaHead velocity', 'to Purkinje cell'},'Interpreter','Tex');
            ylabel(hs(2),{'\DeltaHead velocity', 'to brainstem'},'Interpreter','Tex');
            my_export_fig(fullfile(I.figures_path,sprintf('learnDeltaW_%s_PF%g.pdf', types{jj}, PFs(ii)) ))
        end
    end
end




%% Plot change in vestibular weights after learning - ALL FEEDBACK STRENGTHS
ymax1 = 0.2;
types = {'LOW','HIGH'};

for jj =1:length(types)
    figure('Units','centi','Pos',pos_learn)
    
    
    % Normal learning
    K_baseline = K1;
    if strcmp(types{jj},'LOW');     K_learn = K0;
    else                            K_learn = K2;
    end
    
    for ii = 1:length(PFs)
        
        if length(K_learn)<ii || isempty(K_learn(ii).PH); continue; end
        hs(1) = subplot(2,1,1);
        plot(tt_H(tt_H<=maxt+eps), K_learn(ii).PH(tt_H<=maxt+eps) - K_baseline(ii).PH(tt_H<=maxt),'Color',I.c(ii,:),'Clipping','off'); hold on
        
        hs(2) = subplot(2,1,2);
        plot(tt_H(tt_H<=maxt+eps), -(K_learn(ii).EH(tt_H<=maxt+eps) - K_baseline(ii).EH(tt_H<=maxt)),'Color',I.c(ii,:),'Clipping','off'); hold on
        
    end
    
    axis(hs,'square')
    
    set(hs,'XLim',[0 maxt+dt],'XTick',0:25:maxt, 'Box','off')
    set(hs(1:2), 'YLim',[-ymax1 ymax1],'YTick',-ymax1:dy:ymax1)
    subplot(hs(1))
    h = plot(tt_H,zeros(size(tt_H)),'-','Color',.5*[1 1 1]); uistack(h,'bottom');
    subplot(hs(2))
    h = plot(tt_H,zeros(size(tt_H)),'-','Color',.5*[1 1 1]); uistack(h,'bottom');
    xlabel('Time (ms)')
    
    if option_final
        fixticks(ptick)
        cleanticks;
         export_fig(fullfile(I.figures_path,sprintf('learnDeltaW_%s_all.pdf', types{jj})))
    else   
        title(hs(1), types{jj})
        ylabel(hs(1), {'\DeltaHead velocity', 'to Purkinje cell'},'Interpreter','Tex');
        ylabel(hs(2),{'\DeltaHead velocity', 'to brainstem'},'Interpreter','Tex');
        my_export_fig(fullfile(I.figures_path,sprintf('learnDeltaW_%s_all.pdf', types{jj})))
    end
end



%% Plot change in vestibular weights after learning - NET
ymax1 = 0; ymax2 = 0;
hs = [];
dy = .2;
bar_width = 0.08;

types = {'Low','High'};
for jj = 1:length(types)
    figure('Units','centi','Pos',pos_learn)
    
    % Normal learning
    K_baseline = K1;
    if strcmpi(types{jj},'Low');     K_learn = K0;
    else                             K_learn = K2;
    end
        
    for ii = 1:length(PFs)
        
        if length(K_learn)<ii || isempty(K_learn(ii).PH); continue; end
        
        delta_PH = K_learn(ii).PH - K_baseline(ii).PH;
        delta_EH = K_learn(ii).EH - K_baseline(ii).EH;
        
        hs(1) = subplot(2,1,1);
        bar(I.PFs(ii), sum(delta_PH), bar_width, 'FaceColor',I.c(ii,:),'EdgeColor','none','Clipping','off'); hold on;
        
        hs(2) = subplot(2,1,2);
        h_temp = bar(I.PFs(ii), -sum(delta_EH), bar_width, 'FaceColor',I.c(ii,:),'EdgeColor','none','Clipping','off'); hold on;
        
        ymax1 = max(ymax1, max(ceil(abs(sum(delta_PH))/dy))*dy);
        ymax2 = max(ymax2, max(ceil(abs(sum(delta_EH))/dy))*dy);
        
    end
    
    set(hs, 'YLim',[-ymax1 ymax1],'YTick',-ymax1:dy:ymax1)    
    set(hs, 'XLim',[min(PFs)-bar_width max(PFs)+bar_width],'XTick',sort(PFs)) % If only plotting two        
            xlabel('Feedback')

    if option_final
        axis(hs,'square')
        fixticks(ptick); drawnow;
        cleanticks
        cleanticks('x')
        export_fig(fullfile(I.figures_path,sprintf('learnDeltaW_%s_net.pdf', types{jj})),'-nocrop')
    else
        title(hs(1), types{jj})
        ylabel(hs(1), {'\DeltaHead velocity', 'to Purkinje cell'},'Interpreter','Tex');
        ylabel(hs(2),{'\DeltaHead velocity', 'to brainstem'},'Interpreter','Tex');
        my_export_fig(fullfile(I.figures_path,sprintf('learnDeltaW_%s_net.pdf', types{jj})))
    end
end


%% Replot steps after learning
for ii = [1 11]
    hF_step = figure;   drawnow; 
    plotStep(K1(ii), K2(ii), K0(ii), I, I.learn_color)
    if option_final
        my_export_fig(hF_step, fullfile(I.figures_path, sprintf('learn_step%g.pdf', I.PFs(ii))));
    end
end


% %% Plot cancellation traces
% for ii = [1 11]
% plotCancSplitTrace(I, K1(ii), K2(ii), K0(ii));
% end



%% Revert to factory default
set(groot,'defaultAxesXTickLabelRotationMode','remove')
set(groot,'defaultAxesYTickLabelRotationMode','remove')
set(groot,'defaultAxesZTickLabelRotationMode','remove')


%}