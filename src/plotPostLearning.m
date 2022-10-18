function [hF_freq_eye, hF_freq_pc, hF_filters, hF_step, hF_sine] = plotPostLearning(...
    I, K1, K2, K0, RL_data, tts, head, conds)
% Helper function to create a bunch of plots after running VOR model  learning


% Get frequency response after learning
freqs = RL_data.freqs; % Use frequencies from Ramachandran & Lisberger 2005
% freqs = logspace(log10(0.5), log10(50),25); % Use more continuous frequencies
[~, ~, E_gain1, E_phase1, P_gain1, P_phase1] = getFreq(K1, I, freqs);
[~, ~, E_gain2, E_phase2, P_gain2, P_phase2] = getFreq(K2, I, freqs);
[~, ~, E_gain0, E_phase0, P_gain0, P_phase0] = getFreq(K0, I, freqs);

% Plot eye frequency response after learning
hF_freq_eye = figure; pause(.2); drawnow;
plotFreq(freqs, E_gain1*RL_data.scaleJR2SL, E_phase1, 'o-', I.learn_color(2,:), I.learn_color(2,:));
plotFreq(freqs, E_gain2*RL_data.scaleJR2SL, E_phase2, 'o-', I.learn_color(3,:), I.learn_color(3,:));
plotFreq(freqs, E_gain0*RL_data.scaleJR2SL, E_phase0, 'o-', I.learn_color(1,:), I.learn_color(1,:));
xlabel('Frequency'); subplot(211); title('Eye')

% Plot pc frequency response after learning
hF_freq_pc = figure; pause(.2); drawnow;
plotFreq(freqs, P_gain1, P_phase1, 'o-', I.learn_color(2,:), I.learn_color(2,:));
plotFreq(freqs, P_gain2, P_phase2, 'o-', I.learn_color(3,:), I.learn_color(2,:));
plotFreq(freqs, P_gain0, P_phase0, 'o-', I.learn_color(1,:), I.learn_color(1,:));
xlabel('Frequency'); subplot(211); title('PC')

% Plot filters before and after learning
hp_filter  = [];
hF_filters = figure; hold on; pause(.2); drawnow;
[~, hp_filter([2 3])] = plotFiltersLearn(K1, K2, I.learn_color(2,:), I.learn_color(3,:));
[~, hp_filter([2 1])] = plotFiltersLearn(K1, K0, I.learn_color(2,:), I.learn_color(1,:));
legend(hp_filter, {'Low','Normal','High'},'Box','off','Location','SouthEast')
linkaxes(findobj(gcf,'Type','axes'),'x'); xlim([0 I.max_time_basis])

% Plot step response
hF_step = figure; drawnow; plotStep(K1, K2, K0, I, I.learn_color)
% Plot sine responses
hF_sine = figure; drawnow; plotSine(K1, K2, K0, I, I.learn_color)
 