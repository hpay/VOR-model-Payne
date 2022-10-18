function [hLine, hFill] = fillerror(x, ymean, yerr, color, lineon, face_alpha,  varargin)
% FILLERROR plot filled error bars and optional central line
% [h hFill] = FILLERROR(x,mean,error,c,lineon, lighten,varargin)
%
% INPUTS
% x - vector of x positions to plot
% mean - vector of y positions to plot (usually the mean)
% error - vector containing the width of error bars (usually SEM), or a 2 x
% n matrix containing positive and negative error bars
% c - (optional) color ('k' default). The fill can be lightened.
% lineon - (optional) set to 1 to include the mean line (default 1)
% lighten - (optional) fraction from 0 to 1 to lighten the background
% Include other param value pairs to format the line or use the handle outputs
%
% example:
% x = 1:10;
% ymean = rand(1,10);
% yerr = rand(1,10);
% [hLine, hFill] = fillerror(x, ymean, yerr)
%
% Author: Hannah Payne
% Date: 2012

% Defaults
if ~exist('lineon','var')  || isempty(lineon);  lineon  = 1; end
if ~exist('color','var') || isempty(color);   color = 'k';   end
if ~exist('face_alpha','var') || isempty(face_alpha); face_alpha = 0.5; end

if ~isempty(varargin) && strcmp(get(varargin{1}, 'type'), 'axes')
    ax = varargin{1};
    varargin = varargin(2:end);
else
    ax = gca;
end

% Columnize everything
x = x(:);
ymean = ymean(:);
if size(yerr,1)<size(yerr,2)
    yerr = yerr';
end
if size(yerr,2)==1
    yerr = [-abs(yerr) yerr];
end

% Get rid of any NaNs in the error bars
yerr(isnan(yerr)) = 0;

% Create vectors
plotx = [x; flipud(x)];
ploty = [ymean+yerr(:,1); flipud(ymean+yerr(:,2))];
ploty(isnan(ploty)) = 0;

% Plot it
hFill = patch(ax, plotx,ploty,color,'EdgeColor','none',varargin{:});
set(hFill,'FaceAlpha',face_alpha)

% Plot line
if lineon
    hLine = line(ax, x,ymean, 'Color',color,varargin{:});
else
    hLine = [];
end
