function cout = lighten(cin, factor)
% function cout = lighten(cin, factor)
% Default 50%
% Hannah Payne

if ~exist('factor','var')
    factor = .5;
end

if ischar(cin) 
    if strcmp(cin, 'none')
    cout = 'none';
    return
    else
        cin = colorspec(cin);
    end
end


cout = (1-cin)*factor+cin;
