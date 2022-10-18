function my_export_fig(varargin)
% my_export_fig
% my_export_fig((f_handle), (savename), (type), (pad))
%
% Make matlab print to the correct size of the figure on the screen
%
% pad can be positive or negative to either increase or decrease the
% whitespace around the figure
%
% type is usually unnecessary if extension is provided
%
% Hannah Payne
% 9/2017

% Set defaults
f_handle = gcf;
filename = 'figure.pdf';
type = [];
pad = 0;

for ii = 1 : length(varargin)
    cur_arg = varargin{ii};
    
    if isempty( cur_arg )
        
    elseif ~ischar( cur_arg )
        % Non-string argument better be a handle of a Figure or model
        f_handle = cur_arg;
        
    elseif (cur_arg(1) ~= '-') % Filename
        filename = cur_arg;
        
    else
        switch( cur_arg(2) )
            case 'd'
                type = cur_arg;
                
            case 'p' % Padding
                pad = str2double(cur_arg(3:end));
        end
        
    end
end

set(f_handle,'Units','Inches');
pos = get(f_handle,'Position');
set(f_handle,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize', (1+pad)*[pos(3), pos(4)])

[~, ~, fileterm] = fileparts(filename);

if isempty(type)
    switch  fileterm
        case '.pdf'
            type = '-dpdf';
            option = '-painters';
            
        case {'.jpg' '.jpeg'}
            type  = '-djpeg';
        case '.png'
            type = '-dpng';
        case {'.tif', 'tiff'}
            type = '-dtiff';
        case '.bmp'
            type = '-dbmp';
        case '.eps'
            type = '-deps';
            option = '-painters';
            
        case '.emf'
            type = '-demf';
            option = '-painters';
            
        case '.svg'
            type = '-dsvg';
            option = '-painters';
            
        case '.ps'
            type = '-dpsc';
            option = '-painters';
            
        otherwise
            type = '-dpdf';
            option = '-painters';
    end
end
try
if exist('option','var')
    print(f_handle, filename, type,option)
else
    print(f_handle, filename, type)
end
catch
    warning('Couldn''t print!')
end
