clc
clear all
close all

%%
t = deg2rad(0:10:180);
p = deg2rad(0:10:359);

[t,p] = meshgrid(t,p);

%
x = sin(t).*cos(p);
y = sin(t).*sin(p);
z = cos(t);

n = [1, 1, 1];

%%
plot(x,z, 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'y');
axis equal
set(gca, 'visible', 'off')