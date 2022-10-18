function  h = plotInputsModel(fHandle, times, ipsiEye, ipsiPC,color, offsets, xlims, subplots, varargin)
% Just plot model fits
% NEW: single graph for all conditions
figure(fHandle);

if ~exist('subplots','var')
    subplots(2) = subplot(1, 3, 2);
    subplots(3) = subplot(1, 3, 3);
end

h = [];
for ii = 1:length(ipsiEye)           
        
    % Scale x axis to match
    xlims_actual = get(subplots(2),'XLim');
    mask = times{ii}>=xlims(ii, 1) & times{ii}<=xlims(ii,2);
    ttPlot = scaledata(times{ii}(mask), xlims_actual);  
   
    % Plot PC
    h(ii,1) = line(subplots(2), ttPlot, offsets(ii) + nanmean(ipsiPC{ii}(mask,:),2),'Color',color,'LineWidth',1,'Clipping','off',varargin{:});
    
    % Plot Eye
    h(ii,2) = line(subplots(3), ttPlot,offsets(ii) +nanmean(ipsiEye{ii}(mask,:),2),'Color',color,'LineWidth',1,'Clipping','off',varargin{:});
end
drawnow; % to allow plot to display
