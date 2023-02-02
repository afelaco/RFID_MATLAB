function [G_R_T, G_T_R] = importS2P_2(f, Folder)

selpath = uigetdir(Folder);
addpath(selpath)

listing = dir([selpath '\*.s2p']);

for n = 1 : length(listing)
    
    sobj_R_T = sparameters(listing(n).name);
    sobj_T_R = sparameters(flip(flip(sobj_R_T.Parameters, 1), 2), sobj_R_T.Frequencies, sobj_R_T.Impedance);
    
    [~,I] = min(abs(sobj_R_T.Frequencies-f));
    
    idx = sscanf(listing(n).name, '%d_%d');
    
    i = idx(1);
    j = idx(2);
    
    if i == 0
        
        G_R_T.mid(j) = 20.*log10(abs(sobj_R_T.Parameters(2, 1, I)));
        
        G_T_R.mid(j) = 20.*log10(abs(sobj_T_R.Parameters(2, 1, I)));
        
    else
        
        G_R_T.grid(i,j) = 20.*log10(abs(sobj_R_T.Parameters(2, 1, I)));
        
        G_T_R.grid(i,j) = 20.*log10(abs(sobj_T_R.Parameters(2, 1, I)));
        
    end
    
end

end