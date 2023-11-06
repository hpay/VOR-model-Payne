function [dHandle, sHandle] = plotFreq(frequencies, gains, phases,...
linespec,color,fillcolor, sHandle, newplot)
% Make bode plots with data and model output
% plotFreq(RL_data.freqs, RL_data.gains1, RL_data.phases1+180,
% linespec, I.learn_color(2,:), I.learn_color(2,:));

lw = 1.2;
msize = 3;
phases = unwrap(wrapTo180(phases-180)*pi/180)*180/pi;
linestyle = '-';
marker = 'o';
if ~exist('color','var')
    color = 'k';
    fillcolor = 'k';
end

if exist('linespec','var')
    
    if ~strfind(linespec, '-')
        linestyle = 'none';
    else
        linespec(regexp(linespec,'[-.]'))=[];
        marker = linespec;
        if isempty(marker); marker = 'none'; end
    end
end

if ~exist('newplot','var')
    newplot = 0;
end

if exist('sHandle','var') && ~newplot 
    axes(sHandle(1))
    dHandle(1) = line(log10(frequencies),log10(gains),'Marker', marker,...
        'Color',color,'MarkerEdgeColor','none','MarkerFaceColor',fillcolor,...
        'MarkerSize',msize,'LineStyle',linestyle,'Clipping','off','LineWidth',lw); hold on % DATA
    axes(sHandle(2))
    dHandle(2) = line(log10(frequencies), rad2deg(unwrap(deg2rad(wrapTo180(phases)))),'Marker', marker,...
        'Color',color,'MarkerEdgeColor','none','MarkerFaceColor',fillcolor,...
        'MarkerSize',msize,'LineStyle',linestyle,'Clipping','off','LineWidth',lw); hold on % DATA
else        
    ylimsPhase = [-100 30];
    if ~exist('sHandle','var')
        hold on;
        sHandle = subplot(2,1,1); hold on;
    else
        axes(sHandle(1)); hold on
    end
    htemp = line(log10([frequencies(1) frequencies(end)]),log10([1 1]),'Color','k','LineStyle',':'); hold on; uistack(htemp, 'bottom')
    dHandle(1) = plot(log10(frequencies),log10(gains), 'Marker', marker,'Color',color,...
        'MarkerEdgeColor','none','MarkerFaceColor',fillcolor,'MarkerSize',msize,'LineStyle',linestyle,'LineWidth',lw); hold on % DATA

    set(gca,'YTick',log10([.1:.1:1 2 3]),'YTickLabel',{[],[],[],[],'.5',[],[],[],[],'1',[],'3'})
    ylim(log10([.3 3]));
    set(gca,'XTickLabel','')    
    
    if length(sHandle)<2
        hold on;
        sHandle(2) = subplot(2,1,2); hold on;
    else
        axes(sHandle(2)); hold on
    end
    htemp = line(log10([frequencies(1) frequencies(end)]),[0 0],'Color','k','LineStyle',':'); hold on; uistack(htemp, 'bottom')
    
    dHandle(2) = plot(log10(frequencies),rad2deg(unwrap(deg2rad(wrapTo180(phases)))), 'Marker', marker,...
        'Color',color,'MarkerEdgeColor','none','MarkerFaceColor',fillcolor,'MarkerSize',msize,'LineStyle',linestyle,'LineWidth',lw); hold on % DATA

    ylim(ylimsPhase);
    yticks = -90:30:ylimsPhase(2);
    temp = arrayfun(@(x) [num2str(x) '\circ'], yticks,'UniformOutput',false);
%     temp(ismember(yticks, [-60  -30])) = {''};
    set(gca,'YTick',yticks,'YTickLabel',temp)
    
    set(sHandle,'XLim',log10([frequencies(1) frequencies(end)]),'Box','off')
    
    xticks = [.5:.1:1 2:10 20:10:frequencies(end)];
    temp = arrayfun(@num2str,xticks,'UniformOutput',false);
    xticklabels = cell(size(temp)); inds = [1 6 10 16 length(xticks)];
    xticklabels(inds)=temp(inds);
    set(sHandle,'XTick',log10(xticks))
    set(gca,'XTickLabel',xticklabels)
    
    set(findobj(gcf,'Type','Patch','-or','Type','Line'),'Clipping','off')
end