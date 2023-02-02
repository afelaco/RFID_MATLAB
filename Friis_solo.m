clc
clear all
close all

%% Import transmitter.
TX = getdut('TX');

%% Import receiver.
RX = getdut('RX');

%% Lowest order.
L = min(TX.L, RX.L);
d = max(TX.d, RX.d);
w = TX.w;
k = TX.k;

%%
slater = load('Lib\Slater.mat');
slater.MM = slater.MM(1:(L+1)^2, 1:(L+1)^2, 1:(L+1)^2);
slater.ME = slater.ME(1:(L+1)^2, 1:(L+1)^2, 1:(L+1)^2);

%% Grid.
Rho_min = 0.5*w;
Rho_max = 3*w*sqrt(2);

Grid.Theta = unique(TX.F.Theta);
Grid.Phi = unique(TX.F.Phi);
Grid.Rho = linspace(Rho_min, Rho_max);

[Grid.Theta, Grid.Phi, Grid.Rho] = ndgrid(Grid.Theta, Grid.Phi, Grid.Rho);

%% SSH.
Y = SSH(L, Grid);

%% TX VSHA.
[TX.VSHA.e, TX.VSHA.m] = VSHA(TX.F, L);

% viewCoefficients(e, "$e$")
% viewCoefficients(m, "$m$")

%% RX VSHA.
[RX.VSHA.e, RX.VSHA.m] = VSHA(RX.F, L);

% viewCoefficients(e, "$e$")
% viewCoefficients(m, "$m$")

%% RX VSHA inversion.
i = (0:(L+1)^2-1)';
l = floor(sqrt(i));

[E, M] = VSH(Y);

F_TX = VSHS(TX.VSHA.e, TX.VSHA.m, E, M);

F_RX = VSHS((-1).^(l+1).*RX.VSHA.e, (-1).^(l).*RX.VSHA.m, E, M);

%% Voltage Friis.
Z_L = 50;

RX.I = RX.V./(RX.Z + Z_L);

V_L_Friis = -2.*1i.*w./(120.*pi.*Grid.Rho.*RX.I).*dot(F_TX, F_RX, 3).*(Z_L./(RX.Z + Z_L));