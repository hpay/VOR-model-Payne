function cleanticks(ax, nkeep, ah)
% CLEANTICKS clean up tick labels
%  CLEANTICKS   removes central tick labels from y axis 
%  CLEANTICKS('y')   removes central tick labels from y axis 
%  CLEANTICKS('x') removes central tick labels from x axis 
%  CLEANTICKS('x', nkeep) only keeps every nth tick
%  CLEANTICKS('x', [nkeep nstart]) only keeps every nth tick and starts at
%  nstart

if ~exist('ax','var') || isempty(ax)
    ax = 'y';
end

if strcmpi(ax,'x')
    TickLabels = 'XTickLabel';
    Tick = 'XTick';
    Lim = 'XLim';
elseif strcmpi(ax,'y')
    TickLabels = 'YTickLabel';
    Tick = 'YTick';
    Lim = 'YLim';
else
    strcmpi(ax,'z')
    TickLabels = 'ZTickLabel';
    Tick = 'ZTick';
    Lim = 'ZLim';
end

if ~exist('ah','var')
ah = findobj(gcf,'Type','Axes');
end
for i = 1:numel(ah)
    
    % Make sure ticks correspond to those actually showing
        lims = get(ah(i),Lim);

    origTicks = get(ah(i), Tick);
    newTicks = origTicks;
    newTicks(origTicks<lims(1) | origTicks>lims(2)) = [];
    set(ah(i),Tick,newTicks);
    
    origLabel = get(ah(i), TickLabels);
    if ~iscell(origLabel)
        origLabel = mat2cell(origLabel, ones(size(origLabel,1),1),size(origLabel,2));
    end
    newLabel = origLabel;
    
    if exist('nkeep','var') && length(nkeep)>1
        nstart= nkeep(2);
        nkeep = nkeep(1);
    else        
        nstart= 1;
    end
    
    if ~exist('nkeep','var') || isempty(nkeep)
        newLabel(2:end-1) = {''};
        % Always keep '0'
        newLabel(strcmp(origLabel,'0')) = {'0'};
    else
        newLabel(1:end) = {''};
        newLabel(nstart:nkeep:end) = origLabel(nstart:nkeep:end);
    end    
 
    set(ah(i), TickLabels, newLabel)
    
end