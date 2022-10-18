function shrink(f)
% Shrink window by fraction from  default size
% shrink % defualt 50%
% shrink(f) % shink by f
% shrink([h v])

if ~exist('f','var')
    f = [.5 .5];
end

if length(f)==1;
    f = [f f];
end
tempUnits = get(gcf, 'Units');
set(gcf, 'Units', 'pixels');

a = get(gcf,'Position');
set(gcf,'Position',[a(1) a(2) a(3)*f(1) a(4)*f(2)])
drawnow;
pause(.1)
set(gcf, 'Units', tempUnits);
