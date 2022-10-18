function [tts, head, target,  hevel, hepos_err, PC, sines, lights, dt, n_cells, cond_inds] =...
    loadJR(data_path, conds, t_step, flag_average)
% varargout = loadJR(conds)
% [times, head, target,  ipsiEye, ipsiEye_pos, ipsiPC, sines, lights, dt, n_cells, cond_inds] = loadJR(conds)

if ~exist('flag_average','var')
    flag_average = 1;
end


if ~exist('conds','var') || isempty(conds)
    % ------------ JR DATA FOR MODEL FIT -----------------
    conds = {'05Hz_dark'; '05Hz_pursuit';'05Hz_x2';'05Hz_x0';...
        '2Hz_dark'; '2Hz_x2';'2Hz_x0';...
        '5Hz_dark'; '5Hz_x2';'5Hz_x0';...
        '10Hz_dark'; '10Hz_x2';'10Hz_x0';...
        '80ms_dark';'80ms_x2';'80ms_x0';...
        '150ms_dark';'150ms_x2';'150ms_x0';...
        '250ms_dark';'250ms_x2';'250ms_x0';...
        '500ms_dark';'500ms_x2';'500ms_x0'};
end
    

% Path where data is stored
data_path = fullfile(data_path, 'JR_DATA_100HzEyeFilter');

% time to crop steps to - average of ipsi and contra responses
crop = t_step + 0.05;   % set to 2 s - tBefore to equal sine waves timescale
delay_PC = 0.001;       % (ms) delay for window discriminator to trigger spikes

% DOWNSAMPLE factor (from ,5 to 2 ms)
upsamp = 1;
ps = upsamp; % Integer, increase to upsample
qs = 1; % Integer, increase to downsample

warning('off','stats:regress:NoConst')

n = length(conds);

xlims = [];

% Identify dark conditions and sine freqs
lights = cellfun(@isempty,strfind(conds,'dark'));
sines = double(cellfun(@isempty,strfind(conds,'ms'))); % freq if sine, 0 if step
for ii = 1:n
    temp = conds{ii}(1:strfind(conds{ii},'Hz')-1);
    if isempty(temp); sines(ii) = 0; else temp(temp(1)=='0')='.';
        sines(ii)= str2double(temp); end
end

for ii = 1:n
    
    % Sensory inputs
    A = load(fullfile(data_path, [conds{ii} '_head']));
    ttOld = A.time;
    dtOld = mean(diff(ttOld));
    allcellsipsi = A.allcellsipsi;
    
    dt = dtOld*qs/ps;
    
    ttNew = (ttOld(1):dt:ttOld(end))';
    
    tts{ii} = ttNew;
    
    head{ii} = nanmean(allcellsipsi,2);
    head{ii} = interp1(ttOld, head{ii},ttNew,'spline');
    
    load(fullfile(data_path, [conds{ii} '_target']));
    target{ii} = nanmean(allcellsipsi,2);
    target{ii} = interp1(ttOld, target{ii},ttNew,'spline');
    
    % Replace target with perfect sine if its a sinewave
    if sines(ii)
        [amp, phase] = fitsine(tts{ii},head{ii}, sines(ii), 0);
        head{ii} = amp*sin(2*pi*tts{ii}*sines(ii) + 2*pi*phase/360);
        [amp, phase] = fitsine(tts{ii},target{ii}, sines(ii), 0);
        target{ii} = amp*sin(2*pi*tts{ii}*sines(ii) + 2*pi*phase/360);
    end
    
    if ~lights(ii); target{ii}(:) = NaN; end
    
    %% Eye velocity output
    load(fullfile(data_path, [conds{ii}  '_eye']))
    if ps~=1 || qs~=1
        allcellsipsi = interp1(ttOld(:), allcellsipsi, ttNew);
    end
    
    if sines(ii)
        [amp, ph, const] = fitsine(dt,nanmean(allcellsipsi,2),sines(ii),1,0); % remove any baseline for sines
        allcellsipsi = allcellsipsi - const;
    end
    
    if flag_average
        hevel{ii} = nanmean(allcellsipsi,2);
    else
        hevel{ii} = allcellsipsi;
    end
    
    %% Eye position error
    load(fullfile(data_path, [conds{ii}  '_posErr']))
    
    if ps~=1 || qs~=1
        allcellsipsi = interp1(ttOld(:), allcellsipsi, ttNew);
    end
    
    if sines(ii)
        [amp, ph, const] = fitsine(dt,nanmean(allcellsipsi,2),sines(ii),1,0); % remove any baseline for sines
        allcellsipsi = allcellsipsi - const;
    end
    
    if flag_average
    hepos_err{ii} = nanmean(allcellsipsi,2);
    else
    hepos_err{ii} = allcellsipsi;        
    end
    
    %% PC simple spike activity
    load(fullfile(data_path, [conds{ii}  '_PC']))
    if ps~=1 || qs~=1
        allcellsipsi = interp1(ttOld(:), allcellsipsi,ttNew);
        allcellscontra = interp1(ttOld(:), allcellscontra,ttNew);
    end
    
    n_cells(ii) = size(allcellsipsi,2);
    
    % Adjust PC rates for delay
    delay_ind = round(delay_PC/dt)+1;

    if sines(ii)
        allcellsipsi = allcellsipsi([delay_ind:end 1:delay_ind-1], :);
        allcellscontra = allcellscontra([delay_ind:end 1:delay_ind-1], :);        
    else
        allcellsipsi = allcellsipsi(delay_ind:end,:);
        allcellscontra = allcellscontra(delay_ind:end,:); 
    end            
    
    if ~crop
        % Don't average ipsi and contra
        if flag_average            
            PC{ii} = nanmean(allcellsipsi,2);
        else
            PC{ii} = allcellsipsi;            
        end
    else
            % Average ipsi and contra
       if flag_average
            PC{ii} =  (nanmean(allcellsipsi,2) -  nanmean(allcellscontra,2))/2;
        else
            ncycles = min(size(allcellsipsi,2), size(allcellscontra,2));
            PC{ii} = (allcellsipsi(:,1:ncycles) -  allcellscontra(:,1:ncycles))/2;
        end
        % Only crop if step
        if ~sines(ii)
            mask_step = tts{ii}<t_step;
            tts{ii} = tts{ii}(mask_step);
            head{ii} = head{ii}(mask_step);
            target{ii} = target{ii}(mask_step);
            hevel{ii} = hevel{ii}(mask_step,:);
            hepos_err{ii} = hepos_err{ii}(mask_step,:);                        
            PC{ii} = PC{ii}(mask_step,:);
        end
        
    end
    
    
    
    
    
    
    
 
 
    
    
    
    
    
    
    
    
    
 
end

cond_inds = {[1 length(tts{1})]};
for ii = 2:length(tts)
    cond_inds{ii} = cond_inds{ii-1}(2) + [1 length(tts{ii})];
end

tts = tts(:);
head = head(:);
target = target(:);
hevel = hevel(:);
hepos_err = hepos_err(:);
PC = PC(:);
warning('on','stats:regress:NoConst')
