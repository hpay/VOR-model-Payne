%% Plot PC activity during VORD and cancellation
function [fH_vord, fH_canc, P_gain_vor, P_gain_canc] =...
    plotVorCanc(I,  K1, K2, K0, flag_gain_up)
% Makes one complete figure for VORD, and one figure for VORC
% One for gain up, and one for gain down


%% Plotting parameters
label_font = 6;
lw = 1.5;
fig_pos = [200   100  150 160];
ytext1 = {'Purkinje cell modulation'; 'sp/s per deg/s'};

I.sine_steadystate = 0;

ylims_canc = [-.8 2.8];
ylims_vord = .8 * [-1 1];

if K1.PF == -1
    ylims_canc = [-3 3];
    ylims_vord = 1.6 * [-1 1];
end

tt =  0:K1.dt:(I.T_steadystate - K1.dt);

% Function to negate gain if abs(phase) > 90 deg
flipGain = @(gain, ph) gain * ((abs(wrapTo180(ph))<90)*2-1);
learnStr = {'Pre', 'Post'};

[P_gain_vor, P_gain_canc] = deal(NaN(3, 1));
[P_phase_vor, P_phase_canc] = deal(NaN(3, 1));
[E_gain_vor, E_gain_canc] = deal(NaN(3, 1));
[E_phase_vor, E_phase_canc] = deal(NaN(3, 1));
plot_conds = 1:2;


if flag_gain_up
    colors = I.learn_color(2:3,:);
else
    colors = I.learn_color([2 1],:);
end
for jj = plot_conds
    
    %% Run for each learning condition
    switch jj
        case 1;         K = K1;
        case 2    
            if flag_gain_up;  K = K2;
            else;             K = K0;
            end

    end
    
    %% Run model - VORD
    freq = 0.5;
    
    % Only fit sine waves to the steady state portion of the response
    mask_ss = (length(tt) - round(1/(K.dt*freq)) + 1):length(tt);
    
    head_curr = sin(2*pi*tt*freq);
    target_curr = zeros(size(tt));
    light_on = 0;
    sine_on = freq; % (Hz);
    [Ehat, Phat] = modelClosedloop(K, I, head_curr, target_curr,...
        light_on,  sine_on);
    
    % Fit a sine wave to the PC and eye response
    [gain, P_phase_vor(jj)]  = fitsine(K1.dt, Phat(mask_ss), freq);
    P_gain_vor(jj) = flipGain(gain, P_phase_vor(jj));
    [gain, E_phase_vor(jj)]  = fitsine(K1.dt, Ehat(mask_ss), freq);
    E_gain_vor(jj) = flipGain(gain, E_phase_vor(jj));
    
    
    %% Run model - cancellation
    head_curr = sin(2*pi*tt*freq);
    target_curr = head_curr;
    light_on = 1;
    [Ehat_canc(:, jj), Phat_canc(:,jj,:)] = modelClosedloop(K, I, head_curr, target_curr,...
        light_on,  sine_on);
    
    % Fit a sine wave to each PC and eye response
    [gain, P_phase_canc(jj)]  = fitsine(K1.dt, Phat_canc(mask_ss,jj), freq);
    P_gain_canc(jj) = flipGain(gain, P_phase_canc(jj));        
    [gain, E_phase_canc(jj)]  = fitsine(K1.dt, Ehat_canc(mask_ss,jj), freq);
    E_gain_canc(jj) = flipGain(gain, E_phase_canc(jj));
    
end


%% Plot VORD
xplot = plot_conds;
axis_pos = [0.35    0.1100    0.55    0.75];
fH_vord = figure('Pos',fig_pos);
set(gca,'Position',axis_pos)

for jj = plot_conds
    bar(xplot(jj), P_gain_vor(jj), .8,'EdgeColor','none',...% ,c_edge_curr(jj,:)
        'FaceColor',colors(jj,:),'LineWidth',lw); hold on;
end
set(gca,'YTick',min(ylims_vord):.4:max(ylims_vord))


ylim(ylims_vord)
xlim([min(xplot(:))-.8 max(xplot(:))+.8])
ylabel(ytext1,'FontSize',label_font)
set(gca,'XTick',plot_conds,'XTickLabel',learnStr)
axis square
fixticks;

%% Plot cancellation 
fH_canc = figure('Pos',fig_pos);
set(gca,'Position',axis_pos)

for jj = plot_conds%  low norm high plot order
    bar(xplot(1,jj), mean(P_gain_canc(jj)),.8,'EdgeColor','none',...
        'FaceColor',colors(jj,:),'LineWidth',lw); hold on;
end

ylim(ylims_canc)
xlim([min(xplot(:))-.8 max(xplot(:))+.8])
ylabel(ytext1,'FontSize',label_font)
set(gca,'XTick',plot_conds,'XTickLabel',learnStr)
axis square
fixticks;



