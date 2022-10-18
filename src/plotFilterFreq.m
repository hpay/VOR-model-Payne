
function plotFilterFreq(K, hf, cs)
if ~exist('hf','var')
    hf = figure;
else
    figure(hf); clf;
end
f = [0.5 1  2  3.33 4 5 8  10  12.5   16.6700   20.   25.   32.5000   40   50];

xlims = [0 0.05];


if ~exist('cs','var')
cs = {'b','k','r','c','m'};
end
sHandle = gobjects(4,3); for k = 1:12; sHandle(k) = subplot(3,4,k); end

for ii = 1:length(K)
    % Plot freq response
    fs = 1/K(ii).dt;
    % f = freqs;
    
    K_EP_PH = conv(K(ii).PH,K(ii).EP);
    h_PH = freqz(K(ii).PH, 1, f,fs);
    h_PH_EP = freqz(K_EP_PH, 1, f,fs);
    h_EH = freqz(K(ii).EH, 1, f,fs);
    K_net = K_EP_PH + [K(ii).EH; zeros(length(K(ii).EP)-1,1)];
    h_net = freqz(K_net, 1, f,fs);
    
    if ii==1; newplot = 1; else newplot = 0; end
    plotFreq(f, abs(h_PH), rad2deg(angle(h_PH)),'o',cs{ii},cs{ii},sHandle(1,1:2),newplot);
    title(sHandle(1),'k_{PH}')
    plotFreq(f, abs(h_PH_EP), rad2deg(angle(h_PH_EP)),'o',cs{ii},cs{ii},sHandle(2,1:2),newplot);
    title(sHandle(2),'k_{PH}*k_{EP}')
    plotFreq(f, abs(h_EH), rad2deg(angle(h_EH)),'o',cs{ii},cs{ii},sHandle(3,1:2),newplot);
    title(sHandle(3),'k_{EH}')
    plotFreq(f, abs(h_net), rad2deg(angle(h_net)),'o',cs{ii},cs{ii},sHandle(4,1:2),newplot);
    title(sHandle(4),'k_{EH} + k_{PH}*k_{EP}'); fixticks
    
    set(sHandle(:,1),'YLim',log10([.1 3]))
    set(sHandle(:,2),'YLim',[-360 30],'YTick',-360:90:30,'YTickLabelMode','auto')
    ylabel(sHandle(1,1),'Gain')
    ylabel(sHandle(1,2),'Phase')
    
    % Plot filters
    ylims = [-1 1]*max(abs([K(ii).PH; K_net]));
    line(sHandle(1,3),(1:length(K(ii).PH))*K(ii).dt, K(ii).PH,'Color',cs{ii});
    line(sHandle(2,3),(1:length(conv(K(ii).PH,K(ii).EP)))*K(ii).dt, conv(K(ii).PH,K(ii).EP),'Color',cs{ii});
    line(sHandle(3,3),(1:length(K(ii).EH))*K(ii).dt, K(ii).EH,'Color',cs{ii});
    line(sHandle(4,3),(1:length(K_net))*K(ii).dt, K_net,'Color',cs{ii});
    
    for k = 1:4; h=line(sHandle(k,3),(1:length(K_net))*K(ii).dt, zeros(1,length(K_net)),'Color','k','LineStyle','--');uistack(h,'bottom');
    end
    
    set(sHandle(:,3),'XLim',[0 length(K(ii).PH)*K(ii).dt],'YLim',ylims)
    set(sHandle(:,3),'XLim',xlims); % ,'YLim',ylims)
    linkaxes(sHandle(:,3))
    arrayfun( @(x) xlabel(x,'Time (s)'),sHandle(:,3));
    
end

set(hf,'units','normalized','outerposition',[0 0 0.9 0.9]);

fixticks;
drawnow


