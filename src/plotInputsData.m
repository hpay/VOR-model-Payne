function [fHandle, subplots, offsets,xlims] = plotInputsData(conds, times, ...
    head, target, ipsiEye, ipsiPC, sines, light, labels, plotMeans, subplots)

d_offset = 4; % spacing between subplots
extra_offset = zeros(1,25);
extra_offset([6:13]) = 10;
extra_offset([5 8 11]) = 30;
extra_offset(14) = 40;
extra_offset([22 25]) = 10;
fHandle = gcf;

if ~exist('subplots','var')
    for ii = 1:3
        subplots(ii) = subplot(1,3,ii); hold on;
    end
end

%% Plot - describe step window here
xlims_actual = [-.05 0.6];

offset = 0;
for ii = 1:length(conds)
    
    %% Mask for plotting
    if sines(ii)
        if sines(ii)==10
            xlims1 = [0 2]; % Stimulus
            xlims2 = [0 1]; % Response
        else
            xlims1 = [0 2];
            xlims2 = [0 2];
        end
    else

        xlims1 = xlims_actual;
        xlims2 = xlims_actual; 
    end
    
    mask1 = times{ii}<=xlims1(2);
    mask2 = times{ii}<=xlims2(2);
    
    ttPlot1 = scaledata(times{ii}(mask1), xlims_actual);
    ttPlot2 = scaledata(times{ii}(mask2), xlims_actual);
    xlims(ii,:) = [min(times{ii}(mask2)) max(times{ii}(mask2))];
    
    % Find offset - max of all 3 columns
    minStim = min(min(head{ii}(mask1)), min(target{ii}(mask1)));
    maxStim = max(max(head{ii}(mask1)), max(target{ii}(mask1)));
    minPC = min(mean(ipsiPC{ii}(mask1,:),2));
    maxPC = max(mean(ipsiPC{ii}(mask1,:),2));
    minEye = min(mean(ipsiEye{ii}(mask1,:),2));
    maxEye = max(mean(ipsiEye{ii}(mask1,:),2));
    
    offset = offset - max([maxStim, maxPC, maxEye]);
    
    offset = offset - d_offset;
    offset = offset - extra_offset(ii);
    
    
    offsets(ii)= offset;
    
    
    %%------ STIMULI ------%%
    if any(abs(head{ii})>.005)
        line(subplots(1), ttPlot1, offset + head{ii}(mask1),'Color','k','LineWidth',1,'Clipping','off');
    end
    if light(ii)
        line(subplots(1), ttPlot1, offset + target{ii}(mask1),'LineStyle','-','Color',[.5 .5 .5],'LineWidth',1,'Clipping','off'); % Change to dashed line in illustrator
    end
    
    
    %%------ PURKINJE CELL -----%%


    % Shaded data
    grey = [.6 .6 .6];
    
    if plotMeans
        % Line data
        hLine = plot(subplots(2), ttPlot2, offset + nanmean(ipsiPC{ii}(mask2,:),2), 'Color',grey);
        hold on
        set(hLine,'Clipping','off')
    else
        % Shaded data - ADD
        [hLine, hFill] = fillerror(ttPlot2, offset + nanmean(ipsiPC{ii}(mask2,:),2),...
            sem(ipsiPC{ii}(mask2,:),2),grey,0,1,subplots(2),'LineWidth',.5);%x,mean,error,c,lineon,lighten
         hold on
        set([hLine, hFill],'Clipping','off')
    end
    
    
    %%------ EYE -----%%

    if plotMeans
        % Line data
        hLine = plot(subplots(3), ttPlot2, offset + nanmean(ipsiEye{ii}(mask2,:),2), 'Color',grey);  hold on
        set(hLine,'Clipping','off')
    else
        % Shaded data
        [hLine, hFill] = fillerror(ttPlot2, offset + nanmean(ipsiEye{ii}(mask2,:),2),...
            sem(ipsiEye{ii}(mask2,:),2),grey,0,1,subplots(3),'LineWidth',.5);  hold on
        set([hLine, hFill],'Clipping','off')
    end
      
    
    % Design left-hand labels
    ytext = conds{ii};
    if ytext(1)=='0'; ytext = ['0.' ytext(2:end)]; end
    if any(ytext=='H'); ytext = [ytext(1:find(ytext=='H')-1) ' ' ytext(find(ytext=='H'):end)]; end
    if any(ytext=='m'); ytext = [ytext(1:find(ytext=='m')-1) ' ' ytext(find(ytext=='m'):end)]; end
    ytext = {ytext(1:strfind(ytext,'_')-1), ytext(strfind(ytext,'_')+1:end)};
    set(gca,'FontSize',7)

    % Update offset minima for next time
    offset = offset + min([minStim, minPC, minEye]);
    
    
end
set(subplots,'Box','off')
linkaxes(subplots)
ylim([offset 0])

set(subplots(1),'XLim',xlims_actual);
set(subplots,'XTick',[0 .5])
set(subplots,'FontSize',8)
set(subplots(2:end),'YColor','w','YTick',[],'XColor','w','XTick',[])
