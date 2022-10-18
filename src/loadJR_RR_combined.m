
function [conds, nConds_JR, tts, head, target, hevel, ...
    PC, sines, light, dt, n_cells, RR_data] = ...
    loadJR_RR_combined(I, data_path)

% % ------------ JR DATA FOR MODEL FIT -----------------
conds = {'05Hz_dark'; '05Hz_pursuit';'05Hz_x2';'05Hz_x0';...
    '2Hz_dark'; '2Hz_x2';'2Hz_x0';...
    '5Hz_dark'; '5Hz_x2';'5Hz_x0';...
    '10Hz_dark'; '10Hz_x2';'10Hz_x0';...
    '80ms_dark';'80ms_x2';'80ms_x0';...
    '150ms_dark';'150ms_x2';'150ms_x0';...
    '250ms_dark';'250ms_x2';'250ms_x0';...
    '500ms_dark';'500ms_x2';'500ms_x0'};
nConds_JR = length(conds);

% Load JR's data
[tts, head, target, hevel, heposerr, PC, sines, light, dt, n_cells, cond_inds] = ...
    loadJR(data_path, conds, I.t_step_jr);

% ----------- RAMACHANDRAN & LISBERGER DATA --------
switch I.RR_data_mean_or_W
    case 'mean'
        RR_data = load(fullfile(data_path, 'RamachandranLis2005_DATA_mean'));
 
    case 'meanfix' % Remove the outlier data point (12.5 Hz) to avoid overfitting
        RR_data = load(fullfile(data_path, 'RamachandranLis2005_DATA_mean'));
        mask = RR_data.freqs == 12.5;        
        temp = fieldnames(RR_data);
        for ii = 1:length(temp);     RR_data.(temp{ii})(mask,:,:) = [];   end  
end


% Store gain and phase of JR dark PC sine wave data
ind_temp = find(sines & ~light);
freqs_JR_dark = sines(ind_temp);
freqs_JR = unique(sines(sines>0));
for jj = 1:length(freqs_JR_dark)
    [amp_PC_JR,  phase_PC_JR(jj)] = fitsine(dt, PC{ind_temp(jj)}, sines(ind_temp(jj)));
    [amp_eye_JR,  phase_eye_JR(jj)] = fitsine(dt, hevel{ind_temp(jj)}, sines(ind_temp(jj)));
    [amp_head(jj), phase_head] = fitsine(dt, head{ind_temp(jj)}, sines(ind_temp(jj)));
    gain_PC_JR(jj) = amp_PC_JR/amp_head(jj);
    gain_eye_JR(jj) = amp_eye_JR/amp_head(jj);
    phase_PC_JR(jj) = phase_PC_JR(jj) - phase_head;
    phase_eye_JR(jj) = phase_eye_JR(jj) - phase_head;
end
phase_PC_JR = unwrap(phase_PC_JR*pi/180)*180/pi;
phase_eye_JR = unwrap(phase_eye_JR*pi/180)*180/pi;

% If PC phase is opposite head, make gain negative instead.
if abs(phase_PC_JR(1))>90
    gain_PC_JR = -gain_PC_JR;
    phase_PC_JR = unwrap((phase_PC_JR-180)*pi/180)*180/pi;
end

% Get the sine wave amplitude at 0.5 Hz for JR stimulus (should be 10 deg/s)
sine_amp = round(amp_head(1));

% If the low freq VOR gains across datasets don't match, impossible to fit both
% So define scale factor between JR and SL data: first freq for both is 0.5 Hz
RR_data.scaleJR2SL = RR_data.gains1(1)/gain_eye_JR(1);

% Convert SL gain and phase data to time domain 
freqs = RR_data.freqs;
for jj = 1:length(freqs)
    ind = nConds_JR+jj;
    
    T_temp = round(2*freqs(jj))/freqs(jj); % Max time - should be a multiple of sine period
    tt = dt*(1:round(T_temp/dt))'-dt;
    tts{ind} = tt;
    head{ind} = sine_amp*sin(2*pi*freqs(jj)*tt);
    target{ind} = NaN(size(tt));
    hevel{ind} = -1/RR_data.scaleJR2SL * RR_data.gains1(jj) * sine_amp*sin(2*pi*freqs(jj)*tt + RR_data.phases1(jj)*pi/180);
    
    heposerr{ind} = NaN(size(hevel{ind}));
    
    sines(ind) = freqs(jj);
    light(ind) = 0;
    conds{ind} = sprintf('%gHz_dark',freqs(jj));
    
    % Store the start and end indices for each condition in time
    cond_inds{ind} = cond_inds{ind-1}(2) + [1 length(tt)];
    
    % Interpolate from Jennifers data to create purkinje cell
    % response for intermediate frequencies:
    PC{ind} = NaN(size(tt));
    
    % Check if this frequency is within JR range 0.5 to 10 Hz
    if freqs(jj) <= max(freqs_JR_dark)  && I.data_add_interp_PC
        % If so, estimate using linear interpolation the gain
        % and phase of PC response
        RR_data.gain_PC_SLinterp(jj) = interp1(freqs_JR_dark, gain_PC_JR, freqs(jj));
        RR_data.phase_PC_SLinterp(jj) = interp1(freqs_JR_dark, phase_PC_JR, freqs(jj));
        
        % Generate PC sine wave and store it
        PC{ind} = RR_data.gain_PC_SLinterp(jj) * sine_amp*sin(2*pi*freqs(jj)*tt + RR_data.phase_PC_SLinterp(jj)*pi/180);        
    end
       
end

nConds = length(conds);

