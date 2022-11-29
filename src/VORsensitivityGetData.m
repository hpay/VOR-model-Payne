function d0 = VORsensitivityGetData(filepath, freq)
% For run_sensitivity.m
%
% DAT = VORsensitivityGetData('C:\Dropbox\rlab\model\data\data_100HzEyeFilter', 0.5)
% Gives gain (relative to stimulus) and phase in radians. Negative phase =
% phase lag

if freq == 0.5; freqstr = '05'; else
    freqstr = num2str(freq);
end

types = {'dark','pursuit'};
types_table = {'vord','purs'};
DAT = table;

str = '%sHz_%s_%s.mat';
for ii = 1:length(types)
    if strcmp(types{ii}, 'pursuit')
        R = load(fullfile(filepath, sprintf(str, freqstr, types{ii}, 'target')));
    else
        R = load(fullfile(filepath, sprintf(str, freqstr, types{ii}, 'head')));
    end
    [amp_ref, phase_ref] = fitsine(R.time, R.allcellsipsi, freq);
    
    E = load(fullfile(filepath, sprintf(str, freqstr, types{ii}, 'eye')));
    P = load(fullfile(filepath, sprintf(str, freqstr, types{ii}, 'PC')));
    [amp_E, phase_E] = fitsine(E.time, E.allcellsipsi, freq);
    [amp_P, phase_P] = fitsine(P.time, P.allcellsipsi, freq);
    
    DAT.(types_table{ii})('E_am') = amp_E/amp_ref;
    DAT.(types_table{ii})('E_ph') = deg2rad(phase_E-phase_ref); % check deg or radians
    
    DAT.(types_table{ii})('P_am') = amp_P/amp_ref;
    DAT.(types_table{ii})('P_ph') = deg2rad(phase_P-phase_ref);
     
end

for ii = 1:length(types_table)        
    if abs(DAT.(types_table{ii})('E_ph'))>pi/2
        DAT.(types_table{ii})('E_ph') = wrapToPi(DAT.(types_table{ii})('E_ph')-pi);
        DAT.(types_table{ii})('E_am') = -DAT.(types_table{ii})('E_am');
    end
    if abs(DAT.(types_table{ii})('P_ph'))>pi/2
        DAT.(types_table{ii})('P_ph') = wrapToPi(DAT.(types_table{ii})('P_ph')-pi);
        DAT.(types_table{ii})('P_am') = -DAT.(types_table{ii})('P_am');
    end
end


% Convert to complex notation
d0 = [];
d0.Ed = DAT{'E_am','vord'}*cos(DAT{'E_ph','vord'}) + DAT{'E_am','vord'}*sin(DAT{'E_ph','vord'})*1i;
d0.Pd = DAT{'P_am','vord'}*cos(DAT{'P_ph','vord'}) + DAT{'P_am','vord'}*sin(DAT{'P_ph','vord'})*1i;
d0.Ep = DAT{'E_am','purs'}*cos(DAT{'E_ph','purs'}) + DAT{'E_am','purs'}*sin(DAT{'E_ph','purs'})*1i;
d0.Pp = DAT{'P_am','purs'}*cos(DAT{'P_ph','purs'}) + DAT{'P_am','purs'}*sin(DAT{'P_ph','purs'})*1i;
