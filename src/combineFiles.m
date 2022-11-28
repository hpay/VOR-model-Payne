function combineFiles(filepathname)

[~, savename] = fileparts(filepathname);

PFs = 0:.1:1;

clear E K1 K0 K2
for ii = 1:length(PFs)
    
    try
        D = load(fullfile(filepathname, sprintf('%s_%gPF.mat', savename, PFs(ii))));
    catch
        continue;
    end
    K1(ii) = D.K1;
    K0(ii) = D.K0;
    K2(ii) = D.K2;
    E(ii) = D.E;
    
    % These are the same for all PF strengths:
    I = D.I;
    S = D.S;
    B = D.B;
    R = D.R;
    RR_data = D.RR_data;
    conds = D.conds;
end

save(fullfile(filepathname, sprintf('%s.mat', savename)),'I', 'S', 'B', 'R', 'RR_data', ...
    'K1','K2','K0', 'conds', 'E')
