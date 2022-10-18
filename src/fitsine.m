function [amp, phase, offset, rsquare, residuals, stdresid, trace] =...
    fitsine(time, data, freq, const,ploton)
% [amp, phase (deg), rsquare, residuals, stdresid] = fitsine(time/dt, data, freq, const, ploton)

%% Parameters
if ~exist('const','var') || isempty(const)
    const = 1;
end

if ~exist('ploton','var')
    ploton = 0;
end

%% Fit sine
% Row-ize
time = time(:);
data = data(:);

if length(time)==1
    dt = time;
    time = (1:length(data))'*dt - dt;
end
n = length(data);

keep = ~isnan(data);

y1 = sin(2*pi*freq*time);
y2 = cos(2*pi*freq*time);
constant = ones(n,1); % 1/6/13 Only if you want to allow an offset

if const
    vars = [y1 y2 constant];
else
    vars = [y1 y2];
end
[fit, ~,residuals,~,stats]= regress(data(keep), vars(keep,:));
rsquare=stats(1);
trace = fit'*vars';
trace(~keep) = NaN;
stdresid = std(residuals);

amp = sqrt(fit(1)^2+fit(2)^2);
phase = rad2deg(atan2(fit(2), fit(1)));

if const
    offset = fit(3);
else
    offset = 0;
end

%% Debugging - plot
if  ploton
    clf
    plot(time, data); hold on
    plot(time, trace,'r')
    ylim([min(data) max(data)])
end
