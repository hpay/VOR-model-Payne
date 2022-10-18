function [dataout, B_PE] = filterExponential(tau, dt, datain, sines, delay)

if ~exist('delay','var')
    delay = 0;
end

tt = (0:dt:(tau*5 + delay))';

if tau ==0
    
    % Add a delay
    B_PE = zeros(size(tt));
    B_PE(end) = 1;
    
else
    
    % Filter and add a delay
    B_PE = exp(-(tt-delay)/tau);
    B_PE(tt<delay) = 0;
    B_PE =  B_PE/sum(B_PE);
    
end

dataout = cell(size(datain));
for ii = 1:length(datain)
    nt = length(datain{ii});
    if sines(ii)
        curr_data = [datain{ii}; datain{ii}];
    else
        curr_data = datain{ii};
    end
    
    curr_data_filt = filter(B_PE, 1, curr_data);
    
    if sines(ii)
        curr_data_filt = curr_data_filt(nt+1:end);
    end
    dataout{ii} = curr_data_filt;
    
end
