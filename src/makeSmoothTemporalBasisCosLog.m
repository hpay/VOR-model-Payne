function [bases, param]  = makeSmoothTemporalBasisCosLog(duration, n_bases, dt, ...
    coverage, delay, log_scale, norm_on, plot_on)
%
% Input
%   duration: the time that needs to be covered
%   n_bases: number of basis vectors to fill the duration
%   dt : time step
%   coverage: density of coverage (2 is usually good)
%   plot_on: plot results
%   norm_on: norm peaks
%   delay: min time for the first basis to begin
%   log_scale: Fraction of total duration, decrease to sharpen log drop off
%
% Output
%   bases: structure with basis vectors
%
% Based on JPillow code, modified by Hannah Payne 12/2016
%

% Coverage:  set to 1, 2, 3, etc
if ~exist('coverage','var') || isempty('coverage')
    coverage = 1;
end

if ~exist('delay','var') || isempty('delay')
    delay = 0; % (s)
end

if ~exist('log_scale','var') || isempty('log_scale')
    log_scale= .1;
end

if ~exist('norm_on','var') || isempty('norm_on')
    norm_on = 1;
end


if ~exist('plot_on','var') || isempty('plot_on')
    plot_on = 0;
end

nt = ceil(duration/dt) - 1; % number of bins for the basis functions
ttb = repmat((1:nt)', 1, n_bases); % time indices for basis


%   ^
%  / \
% /   \______
%       ^
%     /   \
% ___/     \___
%              ^
%            /   \
% ________ /       \
% For raised cosine, the spacing between the centers must be 1/4 of the
% width of the cosine


% Code based on http://egocodedinsol.github.io/raised_cosine_basis_functions/

t = 0:dt:duration-dt;
logt = log(t+t(end)*log_scale+eps);

cSt = logt(1);
cEnd = logt( find( t>=(t(end)-delay), 1) );  % Added option to start bases after a delay
db = (cEnd-cSt)/(n_bases+1+coverage);        % Avoid cut-off basis at beginning and ends
c = cSt:db:cEnd;

B = NaN(length(t),n_bases);
for k = 1:n_bases
    k_temp = k+2;
    B(:,k) = (cos(max(-pi, min(pi, 1/coverage*(logt-c(k_temp))*pi/(db) )))+1)/2;
    
    if all(B(:,k)==0)
        B(:,k) = eps;
    end
end

bcenters = exp(c(3:end-2)) - log_scale;

% Account for delay
bcenters = bcenters + delay;
delayi = round(delay / dt);
B = [zeros(delayi, n_bases); B(1:end-delayi,:)];


% Normalize if needed
if norm_on
    for i = 1:n_bases
        B(:,i) = B(:,i)/sum(B(:,i));
    end
end

% Outputs
bases = B;
param.duration = duration;
param.nBases = n_bases;
param.edim = size(B, 2);
param.tr = (ttb - 1)*dt;
param.tcenter = bcenters;

% Plot results
if plot_on
    figure; axis(gca);
    B_plot = B;
    B_plot(B_plot==0) = nan;
    set(gca, 'ColorOrder', gray(size(B,2)+4)); hold on
    plot(param.tr, B_plot);
    xlabel('Time')
    
    hold on
    envelope = nansum(B_plot,2)/coverage;
    plot(param.tr, envelope,'k--');
    ylim([0 max(B(:))])
    
end
