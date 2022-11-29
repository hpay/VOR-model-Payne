function [p, Cost, k5_deltas, PF] = plotCost(fC_matfun, k0, i1, i2, zlims, klabels_plot, complex_on, best_fit_on)
nplot = 151;
apos = [1.5 1.5 .83 .83];

ha = gca;

cla
cs = linspace(-1,1, nplot);
ds = linspace(-1,1, nplot);
[Cs,Ds] = meshgrid(cs,ds);
Cost = NaN(nplot);
opts = optimset('Display','none') ; % For fminsearch.
k5_deltas = zeros(nplot);
PF = zeros(nplot);
for ii = 1:length(cs)
    for jj = 1:length(ds) %
        kd = zeros(size(k0));
        if complex_on
            kd(i1) = cs(ii)*k0(i1)/abs(k0(i1)); % the ratio gives a unit vector in the direction of the original weight
            kd(i2) = ds(jj)*k0(i2)/abs(k0(i2));
            if best_fit_on                
                % find the best kPR that minimizes the function now
                f_fit = @(x) fC_matfun(k0 + kd + [0; 0; 0; 0; x(1)+x(2)*1i]);
                xfit = fminsearch(f_fit, [0 0], opts);
                k5_deltas(jj,ii) = xfit(1)+xfit(2)*1i;
                kd(5) = k5_deltas(jj,ii);
            end
        else
            kd(i1) = cs(ii); % only add on to real component
            kd(i2) = ds(jj);
            if best_fit_on
                f_fit = @(x) fC_matfun(k0 + kd + [0; 0; 0; 0; x]);
                xfit = fminsearch(f_fit, 0, opts);
                k5_deltas(jj,ii) = xfit;
                kd(5) = k5_deltas(jj,ii);
            end
        end
        Cost(jj,ii) = fC_matfun(k0 + kd);
    end
end
% Cost = abs(Cost);
mask1 = Cost<zlims(2);
mask2 = bwareafilt(mask1, 1, 4);
mask = imdilate(mask2, strel('square',3));
Cost(~mask) = NaN;

ax = gca; axis vis3d
colory = Cost;

p = surf(Cs,Ds, Cost,colory, 'EdgeColor','none','FaceColor','interp'); hold on;
xlabel(sprintf('\\Delta%s',klabels_plot{i1}));
ylabel(sprintf('\\Delta%s',klabels_plot{i2}));

grid on
zlabel('Cost')
xlim(ax, [-1 1]); ylim(ax, [-1 1]); zlim(ax, zlims)
axis(ax, 'vis3d')
axis(ax, 'square')
set(ax,'ZScale','linear')

colormap(ax(1),flipud(parula(256)))
set(gca,'CLim',[0 zlims(2)])
lighting gouraud
camlight('headlight')
ha.Position = ha.Position.*apos;
ha.XRuler.TickLabelGapOffset = 0;  % default = +2
ha.YRuler.TickLabelGapOffset = 0;  % default = +2

end

