clc
clear all
close all

%% Y Cut.
addpath('C:\Users\Administrator\OneDrive - UGent\MATLAB\Antennas\Array_NFF\NF\Y Cut')

Xc = 0.078264;
Yc = -0.078264;
Zc = -0.0938524;

for i = 1 : 43
    
    import = readmatrix(sprintf('Y_Cut_step_%d.txt', i));

    NFY.X(44 - i, 1) = unique(import(:,1))./1000 - Xc;
    
    if i == 1
        
        boolean = import(:,3) >= Zc;
        
        NFY.Y = unique(import(:,2))./1000 - Yc;
        NFY.Z = unique(import(boolean,3))./1000;
        
    end

    NFY.Field(44 - i, :, :, 1) = reshape(import(boolean,4) + 1i*import(boolean,5), [105 501]);
    NFY.Field(44 - i, :, :, 2) = reshape(import(boolean,6) + 1i*import(boolean,7), [105 501]);
    NFY.Field(44 - i, :, :, 3) = reshape(import(boolean,8) + 1i*import(boolean,9), [105 501]);
    
end

[NFY.X, NFY.Y, NFY.Z] = ndgrid(NFY.X, NFY.Y, NFY.Z);

%% X Cut.
addpath('C:\Users\Administrator\OneDrive - UGent\MATLAB\Antennas\Array_NFF\NF\X Cut')

for i = 1 : 43
    
    import = readmatrix(sprintf('Y_Cut_step_%d.txt', i));

    NFX.Y(i,1) = unique(import(:,1))./1000 + Yc;
    
    if i == 1
        
        boolean = import(:,3) >= Zc;
        
        NFX.X = unique(import(:,2))./1000 + Xc;
        NFX.Z = unique(import(boolean,3))./1000;
        
    end

    NFX.Field(:, i, :, 1) = reshape(import(boolean,4) + 1i*import(boolean,5), [105 501]);
    NFX.Field(:, i, :, 2) = reshape(import(boolean,6) + 1i*import(boolean,7), [105 501]);
    NFX.Field(:, i, :, 3) = reshape(import(boolean,8) + 1i*import(boolean,9), [105 501]);
    
end

[NFX.X, NFX.Y, NFX.Z] = ndgrid(NFX.X, NFX.Y, NFX.Z);

save('NF','NFX','NFY')