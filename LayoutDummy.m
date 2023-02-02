clc
clear all
close all

W = 30;
H = W;

P = [-W/2, -H/2; W/2, -H/2; W/2, H/2; -W/2, H/2];

DT = delaunay(P);

TR = triangulation(DT, P);

stlwrite(TR,'Layout.stl')