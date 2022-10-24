
function paperDefaults
ill_scale = 1; % 1.28 scaling difference in font sizes in illustrator?
num_font = 6;
label_font = 7;
font = 'Myriad Pro';
set(0,     'DefaultFigureColor', 'white', ...          % Figure background white
    'DefaultAxesColor', 'none', ...
    'DefaultAxesBox', 'off',...
    'DefaultAxesTickDir', 'out',...
    'DefaultTextFontSize', (num_font*ill_scale),...
    'DefaultAxesFontSize', (num_font*ill_scale), ...
    'DefaultColorbarFontSize', (num_font*ill_scale), ...
    'DefaultLegendFontSize', (num_font*ill_scale), ...
    'DefaultAxesFontSizeMode','manual',...
    'DefaultAxesTitleFontSizeMultiplier', label_font/num_font, ... % 10*illscale
    'DefaultAxesLabelFontSizeMultiplier', label_font/num_font, ... % 10*illscale
    'DefaultTextFontName', font,...
    'DefaultAxesFontName', font, ...
    'DefaultAxesTitleFontWeight','normal',...
    'DefaultAxesLineWidth',0.75,...
    'DefaultAxesXColor','k',...
    'DefaultAxesYColor','k',...
    'DefaultFigurePaperPositionMode','auto');