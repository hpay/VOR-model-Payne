function [hf_basefit, rmse_eye, rmse_pc, nrmse_eye, nrmse_pc] = ...
    plotBaselineResults(K, I, conds,  tts, head, target, hevel, PC, light, sines, mask, hf_basefit)
% Plot JR data & model before learning & get error


% Use mask to select conditions if needed
if ~exist('mask','var') || isempty(mask)
    mask = 1:min(I.nConds_JR, length(head));
end

hevel = hevel(mask);
PC = PC(mask);
head = head(mask); 
target = target(mask);
light = light(mask);
sines = sines(mask);
conds = conds(mask);
tts = tts(mask);

% Get model output
[Ehat, Phat] = getModelJR(K, I, head, target, light, sines);

% Calculate error
rmse_eye = cellfun(@(xhat, x) sqrt(mean((xhat - x).^2)), Ehat, hevel);
rmse_pc = cellfun(@(xhat, x) sqrt(mean((xhat - x).^2)), Phat, PC);

% Calculate normalized RMSE
x_range = [range(cell2mat(hevel)) range(cell2mat(PC))];
nrmse_eye =  rmse_eye/x_range(1);
nrmse_pc =  rmse_pc/x_range(2);
fprintf('\nRMS error: eye = %.3f; PC = %.3f\n', nanmean(rmse_eye), nanmean(rmse_pc));
fprintf('NRMS error: eye = %.4f; PC = %.4f\n\n', nanmean(nrmse_eye), nanmean(nrmse_pc));

% Plot data
if ~exist('hf_basefit','var') || isempty(hf_basefit)
hf_basefit = figure;
end

figure(hf_basefit)

% Plot data
[hf_basefit, hs, offsets, xlims] = plotInputsData(conds, tts,...
    head, target, hevel, PC, sines, light,1,1);

% Plot model
plotInputsModel(hf_basefit, tts, Ehat, Phat, ...
    I.c(round(K.PF*10)==round(I.PFs*10),:), offsets, xlims); % Impulse response model
