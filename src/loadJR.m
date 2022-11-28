function [conds, tts, head, target, hevel, PC, sines, light, dt, cond_inds] =...
    loadJR(data_path,  t_step)
% [conds, tts, head, target, hevel, PC, sines, light, dt, cond_inds] = ...
% loadJR(data_path, t_step)
%
% t_step indicates the time at which step data is cropped

% ------------ JR DATA FOR MODEL FIT -----------------
conds = {'05Hz_dark'; '05Hz_pursuit';'05Hz_x2';'05Hz_x0';...
    '2Hz_dark'; '2Hz_x2';'2Hz_x0';...
    '5Hz_dark'; '5Hz_x2';'5Hz_x0';...
    '10Hz_dark'; '10Hz_x2';'10Hz_x0';...
    '80ms_dark';'80ms_x2';'80ms_x0';...
    '150ms_dark';'150ms_x2';'150ms_x0';...
    '250ms_dark';'250ms_x2';'250ms_x0';...
    '500ms_dark';'500ms_x2';'500ms_x0'};

% Path where data is stored
data_path = fullfile(data_path, 'JR_DATA_100HzEyeFilter');

% delay for window discriminator to trigger spikes
delay_PC = 0.001;       % (ms) 

n = length(conds);

% Identify dark conditions and sine freqs
light = ~contains(conds,'dark');
sines = zeros(size(conds)); % Get the frequency of the sine wave
for ii = 1:n
    temp = conds{ii}(1:strfind(conds{ii},'Hz')-1);
    if isempty(temp); sines(ii) = 0; else temp(temp(1)=='0')='.';
        sines(ii)= str2double(temp); end
end

tts = cell(n,1);
head = cell(n,1);
target = cell(n,1);
hevel = cell(n,1);
PC = cell(n,1);

for ii = 1:n
    
    % Sensory inputs
    A = load(fullfile(data_path, [conds{ii} '_head']));
    head{ii} = mean(A.allcellsipsi,2,'omitnan');
    tts{ii} = A.time(:);
    dt = mean(diff(tts{ii}));   
    
    A = load(fullfile(data_path, [conds{ii} '_target']));
    target{ii} = mean(A.allcellsipsi,2,'omitnan');
    
    % Replace head and target with perfect sinusoid to reduce noise
    if sines(ii)
        [amp, phase] = fitsine(tts{ii},head{ii}, sines(ii), 0);
        head{ii} = amp*sin(2*pi*tts{ii}*sines(ii) + 2*pi*phase/360);
        [amp, phase] = fitsine(tts{ii},target{ii}, sines(ii), 0);
        target{ii} = amp*sin(2*pi*tts{ii}*sines(ii) + 2*pi*phase/360);
    end
    
    if ~light(ii); target{ii}(:) = NaN; end
    
    %% Eye velocity output
    A = load(fullfile(data_path, [conds{ii}  '_eye']));
    hevel{ii} = A.allcellsipsi;    
    
    % Remove baseline for sines
    if sines(ii)
        [amp, ph, const] = fitsine(dt,hevel{ii},sines(ii),1,0); 
        hevel{ii} = hevel{ii} - const;
    end            
   
    %% PC simple spike activity
    A = load(fullfile(data_path, [conds{ii}  '_PC']));
    PC{ii} =  (A.allcellsipsi -  A.allcellscontra)/2;  % Average ipsi and contra
    
    % Adjust PC rates for delay
    delay_ind = round(delay_PC/dt)+1;
    if sines(ii)
        PC{ii} = PC{ii}([delay_ind:end 1:delay_ind-1], :);
    else
        PC{ii} = PC{ii}(delay_ind:end,:);
    end            
    
    %% Crop if step
    if ~sines(ii)
        mask_step = tts{ii}<t_step;
        tts{ii} = tts{ii}(mask_step);
        head{ii} = head{ii}(mask_step);
        target{ii} = target{ii}(mask_step);
        hevel{ii} = hevel{ii}(mask_step,:);
        PC{ii} = PC{ii}(mask_step,:);
    end    
 
end

% Time index if all conditions are concatenated together
cond_inds = {[1 length(tts{1})]};
for ii = 2:length(tts)
    cond_inds{ii} = cond_inds{ii-1}(2) + [1 length(tts{ii})];
end

