% Analyze sensitivity to parameters
%
% Hannah Payne 2020
%
% Abbreviations for conditions:
% d = VORD, p = smooth pursuit, s = stim
%
% Abbreviations for inputs:
% H = Head, T = Target, L = Light (1 or 0), E = Eye, P = Purkinje, Z = stim

% PARAMS
stim_on = 0;        % Include stim condition
complex_on = 0;     % Include phase components of data

% Constants
code_root = fileparts(fileparts(which(mfilename)));
addpath(fullfile(code_root, 'src'))
syms H T L E P Z            % Inputs and outputs. (L is Kronecker delta: 1 in light and 0 in dark)
syms Ed Pd Ep Pp Es Ps      % d = VORD, p = smooth pursuit, s = stim
syms k1 k2 k3 k4 k5;        % Unknown parameters (see klabels below)
ks = [k1; k2; k3; k4; k5];
knames = {'k1','k2','k3','k4','k5'};
klabels_plot = {'k_{EH}(ω)', 'k_{EP}(ω)','k_{PH}(ω)','k_{PE}(ω)','k_{PR}(ω)'}';
klabels = {'kEH', 'kEP','kPH','kPE','kPR'}';

% Main equation for eye
fE = E == -k1*H + k2*P;

% Main equation for PC
fP = P ==  k3*H + k4*E + k5*L*(T-H-E) + Z*stim_on;

% Generate equations, one for eye and one for PC in each condition,
% assuming stimuli are known and have unity gain and zero phase.
f_Ed = subs(fE, [H T L Z E P], [1 0 0 0 Ed Pd]);
f_Pd = subs(fP, [H T L Z E P], [1 0 0 0 Ed Pd]);
f_Ep = subs(fE, [H T L Z E P], [0 1 1 0 Ep Pp]);
f_Pp = subs(fP, [H T L Z E P], [0 1 1 0 Ep Pp]);

% Solve equations
if stim_on
    % If exact Ps is specified, model is overfit. So set Z = 1 and solve for Ps  
    Ps = 1./(1-k2*k4); 
    f_Es = subs(fE, [H T L Z E P], [0 0 0 1 Es Ps]);
    f_Ps = subs(fP, [H T L Z E P], [0 0 0 1 Es Ps]);
    S0 = solve([f_Ed, f_Pd, f_Ep, f_Pp, f_Es, f_Ps], ks, 'ReturnConditions',true,'MaxDegree',3);
else
    S0 = solve([f_Ed, f_Pd, f_Ep, f_Pp], ks, 'ReturnConditions',true,'MaxDegree',3);
end

% Get the non-trivial solution (for stim)
S = [];
S.k1 = S0.k1(end);
S.k2 = S0.k2(end);
S.k3 = S0.k3(end);
S.k4 = S0.k4(end);
S.k5 = S0.k5(end);
S.parameters = S0.parameters;
S.conditions = S0.conditions(end);

% Display solution
fprintf('\n--------------------------\nAnalytic solution:\n')
for ii= 1:length(ks)
    fprintf('    %s: %s = %s\n',knames{ii}, klabels{ii},S.(knames{ii}))
end
fprintf('Conditions: \n    '); disp(S.conditions); fprintf('\n')


%% TEST MODEL FIT
fprintf('\n--------------------------\nTest with specific data and kPE:\n\n')

% Load data: gain ("_am") and phase ("_ph") (phase in radians, negative = lag)
d0 = VORsensitivityGetData(fullfile(code_root,'\data\JR_DATA_100HzEyeFilter'), 0.5);
% DAT = VORsensitivityGetData('D:\hannah\Dropbox\rlab\model\data\data_100HzEyeFilter', 0.5);

complexToDC = @(x) abs(x)*(2*(abs(angle(x))<pi/2)-1);
if ~complex_on
  d0.Ed = complexToDC(d0.Ed);
  d0.Pd = complexToDC(d0.Pd);
  d0.Ep = complexToDC(d0.Ep);
  d0.Pp = complexToDC(d0.Pp);
end

% Specify kPE to show one solution
default_kPE = 0.2; 
if stim_on % If stim on, generate expected E response
    d0.Es = double(vpa(subs(solve(default_kPE==S.k4, Es), [Ep Pp], [d0.Ep d0.Pp])));    
end
disp('Data in: ')
disp(d0)
fprintf('kPE = %.2f\n', default_kPE)

% Solve for remaining parameters using the solution found above
fprintf('\nSolved paramters out:\n')
k_struc = [];
if stim_on
    k_struc.EH = subs(S.k1, [Ed Pd Ep Pp Es], [d0.Ed d0.Pd d0.Ep d0.Pp d0.Es]);
    k_struc.EP = subs(S.k2, [Ed Pd Ep Pp Es], [d0.Ed d0.Pd d0.Ep d0.Pp d0.Es]);
    k_struc.PH = subs(S.k3, [Ed Pd Ep Pp Es], [d0.Ed d0.Pd d0.Ep d0.Pp d0.Es]);
    k_struc.PE = subs(S.k4, [Ed Pd Ep Pp Es], [d0.Ed d0.Pd d0.Ep d0.Pp d0.Es]);
    k_struc.PR = subs(S.k5, [Ed Pd Ep Pp Es], [d0.Ed d0.Pd d0.Ep d0.Pp d0.Es]);
else
    k_struc.PE = default_kPE;
    k_struc.EH = subs(S.k1, [Ed Pd Ep Pp], [d0.Ed d0.Pd d0.Ep d0.Pp]);
    k_struc.EP = subs(S.k2, [Ed Pd Ep Pp], [d0.Ed d0.Pd d0.Ep d0.Pp]);
    k_struc.PR = solve(subs(S.k4, [Ed Pd Ep Pp], [d0.Ed d0.Pd d0.Ep d0.Pp]) == k_struc.PE);
    k_struc.PH = subs(S.k3, [Ed Pd Ep Pp S.parameters], [d0.Ed d0.Pd d0.Ep d0.Pp k_struc.PR]);
end
k_vector = [k_struc.EH; k_struc.EP; k_struc.PH; k_struc.PE; k_struc.PR]; % concatenate in vector
k0 = double(vpa(k_vector));  % convert from symbolic to double
for ii = 1:5
    fprintf('%s: %s = %.2f + %.2fi\n', knames{ii}, klabels{ii}, real(k0(ii)), imag(k0(ii)));
end

%% Calculate data using rearranged model eqs so that eye and pc are only on the left-hand-side
fprintf('\nTest -- calculated data:\n')
d1 = [];
[fE_lhs, fP_lhs] = solve([fE, fP], E, P);

% Get hessians and cost function landscape around solution
l = .2; % lighten colors
colors = [.3*[1 1 1]; .6*[1 1 1]; lighten([51 56 142]/255, l);  lighten([245 147 39]/255, l); lighten([189 38 50]/255,l)];
C = (real(E - fE_lhs).^2 + imag(E - fE_lhs).^2) + ...
    (real(P - fP_lhs).^2 + imag(P - fP_lhs).^2);

% Cost function for specific behavioral conditions
if stim_on
    fC = (...
        subs(C, [H T L Z E P], [1 0 0 0 d0.Ed d0.Pd]) + ...     % VORD
        subs(C, [H T L Z E P], [0 1 1 0 d0.Ep d0.Pp]) + ...     % Pursuit
        subs(C, [H T L Z E P], [0 0 0 1 d0.Es Ps]))/3;          % Stim 
else
    fC = (...
        subs(C, [H T L E P], [1 0 0 d0.Ed d0.Pd]) + ...         % VORD
        subs(C, [H T L E P], [0 1 1 d0.Ep d0.Pp]))/2;           % Pursuit
end
gradC = jacobian(fC,ks);
hessC = jacobian(gradC,ks);
fC_matfun = matlabFunction(fC,'vars',{ks});
gradC_matfun = matlabFunction(gradC,'vars',{ks});
hessC_matfun = matlabFunction(hessC,'vars',{ks});
fC_sub = subs(fC, ks, k0);
errFinal = fC_matfun(k0)
gradFinal = gradC_matfun(k0)
hessFinal = hessC_matfun(k0)

%% Plot cost function - brainstem parameters
logscale = 1;
zlims = [0 .25];    % Cost
fpos = [17 15 6 6]; % 3.5 for figures

% Plot cost function for brainstem parameters
fi = figure('Units','centi', 'Pos',fpos);
ax = gca; 
i1 = 1; % first param to plot
i2 = 2; % second param to plot
axes(ax(1)); cla
best_fit_on = 0;
p = plotCost(fC_matfun, k0, i1, i2, zlims, klabels_plot, complex_on, best_fit_on);
title('Brainstem pathways','FontWeight','normal')
views = [20, 18];
set(gca,'XDir','reverse')
view(views(1)+90, views(2))
delete(findobj(gca,'Type','Light'))
camlight(82, -90,'infinite')
camlight('headlight')
set(gca,'ZTick',0:.1:.2)
x = get(gca,'XLabel');
set(x,'Units','norm','Pos',[.95 -.05 0])
y = get(gca,'YLabel');
set(y,'Units','norm','Pos',[.4 -.1 0])
% print('sens_kEH_kEP','-dpng','-r300')

%% Plot cost function - Purkinje cell parameters
fi(2) = figure('Units','centi', 'Pos',fpos);
ax(2) = gca;
i1 = 3;
i2 = 4;
axes(ax(2));
best_fit_on = 1; % Find the kPR that minimizes cost for each combination of inputs
[p, C2, k5_deltas] = plotCost(fC_matfun, k0, i1, i2, zlims, klabels_plot, complex_on, best_fit_on);
title('Purkinje cell pathways','FontWeight','normal')
view(views(1), views(2))
delete(findobj(gca,'Type','Light'))
if stim_on
    camlight(1, -20,'infinite')
else
    camlight(1, 15,'infinite')
end
camlight
drawnow;
set(gca,'ZTick',0:.1:.2)
x = get(gca,'XLabel');
set(x,'Units','norm','Pos',[.4 -.1 0])
y = get(gca,'YLabel');
set(y,'Units','norm','Pos',[.95 -.05 0])
% print('sens_kPH_kPE','-dpng','-r300')   

%% Plot cost function - Purkinje cell parameters kPH and kPR
fi(3) = figure('Units','centi', 'Pos',fpos);
ax(3) = gca;
i1 = 3;
i2 = 5;
axes(ax(3));
best_fit_on = 2; % Find the kPR that minimizes cost for each combination of inputs
[p, C2] = plotCost(fC_matfun, k0, i1, i2, zlims, klabels_plot, complex_on, best_fit_on);
title('Purkinje cell pathways','FontWeight','normal')
view(views(1)+180, views(2))
delete(findobj(gca,'Type','Light'))
if stim_on
    camlight(1, -20,'infinite')
else
    camlight(1, 15,'infinite')
end
camlight
drawnow;
set(gca,'ZTick',0:.1:.2)
x = get(gca,'XLabel');
set(x,'Units','norm','Pos',[.4 -.1 0])
y = get(gca,'YLabel');
set(y,'Units','norm','Pos',[.95 -.05 0])
    
% print('sens_kPH_kPR','-dpng','-r300')