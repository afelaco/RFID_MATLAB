clc
clear all
close all

%% Import data.
Frequency = 3e+09;
Wavelength = physconst('LightSpeed')/Frequency;
k = 2*pi/Wavelength;
L = 7;

%% Grid.
N = 100;
Rho_m = 0.5/sqrt(2);

Grid.Theta = linspace(0, pi, N);
Grid.Phi = linspace(0, 2*pi, 2*N);
Grid.Rho = linspace(0, sqrt(2)*Rho_m, N);

[Grid.Theta, Grid.Phi, Grid.Rho] = ndgrid(Grid.Theta, Grid.Phi, Grid.Rho);

%% SSH.
Y = SSH(L, Grid);

%% SSW.
H = HSW(k, Grid, Y);
J = BSW(k, Grid, Y);

%% Translation.
t = 0;
p = 0;
s = 0.15;

TH = TranslatorH(s, t, p, k, L);
TJ = TranslatorJ(s, t, p, k, L);

n = 1;

% viewField(J{1}, Grid)

c_J = reshape(TJ(n,:), 1, 1, 1, []);

J_T = sum(c_J.*cat(4, J{:}), 4);

viewField(J_T, Grid)