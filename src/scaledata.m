% scaledata.m
% Hannah Payne
% 12/04/2011
%
% SCALEDATA
%   dataout = SCALEDATA(datain) scales datain to the interval [0 1]
%   dataout = SCALEDATA(datain, [minval, maxval]) scales to [minval maxval]
%   dataout = SCALEDATA(datain, minval, maxval) scales to [minval maxval]
%   dataout = SCALEDATA(datain, x1, x2, y1, y2) scales x1 to y1 and x2 to y2 linearly

function dataout = scaledata(datain, x1, x2, y1, y2)

    if nargin == 1
        minval = 0;
        maxval = 1;
    elseif nargin ==2
        minval = x1(1);
        maxval = x1(2);
    elseif nargin ==3
        minval = x1;
        maxval = x2;
    end

    if nargin == 5
        m = (y2-y1)/(x2-x1);
       dataout = m*(datain - x1) + y1;
    else
        dataout = datain - min(datain(:));
        dataout = (dataout/range(dataout(:)))*(maxval-minval);
        dataout = dataout + minval;
    end


    
%         dataout = (datain - min(datain(:)))/range(datain(:));
      

