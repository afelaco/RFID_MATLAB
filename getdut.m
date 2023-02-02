function [output, tag] = getdut(dir)

% Get directory.
selpath = uigetdir(dir);
addpath(selpath)
tag = strsplit(selpath, '\');
tag = tag{end};

% Import frequency.
f = readmatrix('frequency.txt');

% Import voltage.
V = readmatrix('voltage.txt');

% Import impedance.
Z = readmatrix('impedance_mat.txt');
Z = Z(1) + 1i.*Z(2);

% Import layout.
Layout = stlread('mesh.stl');
w = physconst('LightSpeed')/f;
k = 2*pi/w;
d = max(vecnorm(Layout.Points, 2, 2));
L = ceil(k*d);

% Import Far Field.
import = readmatrix('far_field.txt');

[FarField.Phi, FarField.Theta] = meshgrid(deg2rad(unique(import(:,2))), deg2rad(unique(import(:,1))));

FarField.Field(:,:,1) = reshape(import(:,3).*exp(1i.*deg2rad(import(:,4))), size(FarField.Theta))./V;
FarField.Field(:,:,2) = reshape(import(:,5).*exp(1i.*deg2rad(import(:,6))), size(FarField.Theta))./V;

FarField.Theta(:,end+1) = FarField.Theta(:,1);
FarField.Phi(:,end+1) = 2.*pi;
FarField.Field(:,end+1,:) = FarField.Field(:,1,:);

output = struct('f', f, 'V', V, 'Z', Z, 'Layout', Layout, 'w', w, 'k', k, 'd', d, 'L', L, 'F', FarField);

end