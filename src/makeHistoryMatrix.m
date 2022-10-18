function S = makeHistoryMatrix(data, nt_filter, pad_data_or_zero, filtLP)
% function S = make_history_matrix(data, nt_filter, pad_data_or_zero, filtLP)
% Convert vector to matrix for for GLM use
% Hannah Payne
% 
% To do lowpass filter:
% filtLP.dt     % (s)
% filtLP.f      % (Hz)

% Columize
data = data(:);

% Pad either with the last bit of data (if sines) or zeros (steps etc)
% NOTE: before 11/9 was padding with an extra zero, causing the stimulus to
% be shifted 1 time point later than it should
data_pad = [pad_data_or_zero*data(end-nt_filter+2:end); data];

% Filter if needed (use padding)
if exist('filtLP','var') && ~isempty(filtLP)
%     data_pad = smooth(data_pad,filtLP);    
    f_nyquist = (1/(filtLP.dt*2));
    fn = filtLP.f/f_nyquist;
    [b,a] = butter(5,fn);
    
    data_pad= filtfilt(b,a,data_pad);    
end

S = NaN(length(data),nt_filter);
for i = 1:length(data)    
    % Flip so most recent point is to left/first column in matrix
    S(i,:) = flip(data_pad(i:i+nt_filter-1));
end