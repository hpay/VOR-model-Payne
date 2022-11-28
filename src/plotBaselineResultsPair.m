function hf_basefit = plotBaselineResultsPair(K, I, conds, tts, head, target, hevel, PC, light, sines)
% Plot JR data & model before learning & get error
% Use with a pair of inputs K e.g. K1([1 11])

% t_step = 0.7;
% flag_average = 0;
% [tts, head, target,  hevel, ~, PC, sines, lights, dt, n_cells] =...
%     loadJR(I.data_path, [], t_step, flag_average);

I.impulse_or_closedloop = 0; % Check with closed loop

% Plot data
hf_basefit = figure;
xpos_offset = [0 .1 .1 .2 .2];
for ii = 1:5
    subplots(ii) = subplot(1,5,ii); hold on;
    
    % Scale here
    scalex = 1.2;    
    scaley = 1.1;
    pos = get(subplots(ii),'Pos');
    set(subplots(ii),'Pos',[pos(1)+pos(3)*xpos_offset(ii)-pos(3)*(scalex-1)/2 , pos(2)-pos(4)*(scaley-1)/2, pos(3)*scalex, pos(4)*scaley])
    
    
end

for gg = 1:2
    
    % Get model output
    [Ehat, Phat] = getModelJR(K(gg), I, head, target, light, sines);
    
    % Plot data
    [hf_basefit, hs, offsets, xlims] = plotInputsData(conds, tts,...
        head, target, hevel, PC, sines, light,1,0, subplots([1 (0:1)+gg*2]));
    
    % Plot model
    plotInputsModel(hf_basefit, tts, Ehat, Phat, ...
        I.c(round(K(gg).PF*10)==round(I.PFs*10),:), offsets, xlims, subplots([1 (0:1)+gg*2]));
    
    for ii = 1:5
    xlim(subplots(ii), get(subplots(ii),'XLim')+[-.01 0])
    
    end
    
end
