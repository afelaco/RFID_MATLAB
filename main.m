clc
clear all
close all

%% Import transmitter.
[TX, ~] = getdut('C:\Users\Administrator\OneDrive - UGent\RFID\data\tx');

%% Import receiver.
[RX, filetag] = getdut('C:\Users\Administrator\OneDrive - UGent\RFID\data\rx');

%% Lowest order.
L = 10;
d = max(TX.d, RX.d);
w = TX.w;
k = TX.k;

%%
slater = load('Lib\Slater.mat');
slater.MM = slater.MM(1:(L+1)^2, 1:(L+1)^2, 1:(L+1)^2);
slater.ME = slater.ME(1:(L+1)^2, 1:(L+1)^2, 1:(L+1)^2);

%% TX VSHA.
Grid.Theta = TX.F.Theta;
Grid.Phi = TX.F.Phi;
Y = SSH(L, Grid);
[E, M] = VSH(Y);
[TX.VSHA.e, TX.VSHA.m] = VSHA(TX.F, 10);

FieldSph(:,:,1) = exp(-1i.*Grid.Phi)./sqrt(2).*(-cos(Grid.Theta).*TX.F.Field(:,:,1) + 1i.*TX.F.Field(:,:,2));
FieldSph(:,:,2) = exp(1i.*Grid.Phi)./sqrt(2).*(cos(Grid.Theta).*TX.F.Field(:,:,1) + 1i.*TX.F.Field(:,:,2));
FieldSph(:,:,3) = -sin(Grid.Theta).*TX.F.Field(:,:,1);

accuracy = max(mag2db(abs(FieldSph - VSHS(TX.VSHA.e, TX.VSHA.m, E,  M))), [], 'all');

% viewCoefficients(e, "$e$")
% viewCoefficients(m, "$m$")

%% RX VSHA.
[RX.VSHA.e, RX.VSHA.m] = VSHA(RX.F, 10);

% viewCoefficients(e, "$e$")
% viewCoefficients(m, "$m$")

%% RX VSHA inversion.
% i = (0:(L+1)^2-1)';
% l = floor(sqrt(i));
% m = i-l.*(l+1);

% [E, M] = VSH(Y);

% F_TX = VSHS(TX.V.*TX.VSHA.e, TX.V.*TX.VSHA.m, E, M);

% F_RX = VSHS((-1).^(l+1).*RX.V.*RX.VSHA.e, (-1).^(l).*RX.V.*RX.VSHA.m, E, M);

%% Added transmitter data.
Z_L = 50;

% P_a_t = 10^((-45)/10);
% M_TX = 4*real(TX.Z)*real(Z_L)/abs(TX.Z+Z_L)^2;
% U_TX = vecnorm(F_TX, 2, 3).^2./(2*120*pi);
% G_TX = 4*pi.*U_TX./P_a_t;

%% Polarization factor.
% Q = abs(dot(F_TX, F_RX, 3)).^2./((vecnorm(F_TX, 2, 3).^2).*(vecnorm(F_RX, 2, 3).^2));

%% Added receiver data.
% M_RX = 4*real(RX.Z)*real(Z_L)/abs(RX.Z+Z_L)^2;

% P_r = 0.01;
% U_RX = vecnorm(F_RX, 2, 3).^2./(2*120*pi);
% G_RX = 4*pi.*U_RX./P_r;
% rcs = TX.w.^2./(4.*pi).*G_RX.*M_RX.*Q;

%% RX Rotation.
% a = pi/2;
% b = pi/2;
% g = -pi/2;
% 
% if any([a, b, g]) ~= 0
%     
%     D = WignerD(a, b, g, RX.L);
%     
%     RX.VSHA.e = D*RX.VSHA.e;
%     RX.VSHA.m = D*RX.VSHA.m;
%     
% end

%% Grid.
N = 50;
Rho_min = 0;
Rho_max = 3*w*sqrt(2);

Grid.Theta = pi/2;
Grid.Phi = linspace(0, 2*pi, 2*N);
Grid.Rho = linspace(Rho_min, Rho_max, 5*N);

[Grid.Theta, Grid.Phi, Grid.Rho] = ndgrid(Grid.Theta, Grid.Phi, Grid.Rho);

%% SSH.
Y = SSH(L, Grid);

%% SSW.
H = HSW(k, Grid, Y);
% % H = RSW(k, Grid, Y);

%% Downlink.
C_R_T = sum(((TX.VSHA.m)'.*reshape(RX.VSHA.m, 1, 1, []) - ...
    (TX.VSHA.e)'.*reshape(RX.VSHA.e, 1, 1, [])).*slater.MM, [2 3]) + ...
    sum(((TX.VSHA.e)'.*reshape(RX.VSHA.m, 1, 1, []) - ...
    (TX.VSHA.m)'.*reshape(RX.VSHA.e, 1, 1, [])).*slater.ME, [2 3]);

S_R_T = (-1./30).*RX.Z.*TX.V.*sum(reshape(C_R_T, 1, 1, 1, []).*cat(4, H{:}), 4);

V_R_T = S_R_T.*(Z_L./(RX.Z + Z_L));

P_t = 10.*log10(abs(V_R_T).^2./(2.*Z_L)) + 30;

%% Uplink.
% C_T_R = sum(((RX.VSHA.m)'.*reshape(TX.VSHA.m, 1, 1, []) - ...
%     (RX.VSHA.e)'.*reshape(TX.VSHA.e, 1, 1, [])).*slater.MM, [2 3]) + ...
%     sum(((RX.VSHA.e)'.*reshape(TX.VSHA.m, 1, 1, []) - ...
%     (RX.VSHA.m)'.*reshape(TX.VSHA.e, 1, 1, [])).*slater.ME, [2 3]);
% 
% S_T_R = (-1./30).*TX.Z.*S_R_T.*(RX.Z./(RX.Z + Z_L)).*sum(reshape((-1).^l.*C_T_R, 1, 1, 1, []).*conj(cat(4, H{:})), 4);
% 
% V_T_R = S_T_R.*(Z_L./(TX.Z + Z_L));
% 
% RSSI = 10.*log10(abs(V_T_R).^2./(2.*Z_L)) + 30;

%% Downlink Radar.
% L_rt = (4.*pi.*Grid.Rho./TX.w).^2;
% P_t_Radar = 10.*log10(P_a_t.*G_TX.*Q.*M_RX.*G_RX./L_rt) + 30;

%% Uplink Radar.
% S = 4.*pi.*Grid.Rho.^2;
% RSSI_Radar =  10.*log10(P_a_t.*M_TX.*G_TX.^2.*Q.*rcs./(S.*L_rt)) + 30;

%% Interpolate data.
[X, Y, Z] = sph2cart_c(Grid.Rho, Grid.Theta, Grid.Phi);

X = squeeze(X);
Y = squeeze(Y);
P_t = squeeze(P_t);
% RSSI = squeeze(RSSI);
% P_t_Radar = squeeze(P_t_Radar);
% RSSI_Radar = squeeze(RSSI_Radar);

X_q_mid = linspace(1, 3, 5*N).*w;
Y_q_mid = 0;

[X_q_mid, Y_q_mid] = meshgrid(X_q_mid, Y_q_mid);

P_t_q.mid = griddata(X, Y, P_t, X_q_mid, Y_q_mid);
% RSSI_q.mid = griddata(X, Y, RSSI, X_q_mid, Y_q_mid);
% P_t_q_Radar.mid = griddata(X, Y, P_t_Radar, X_q_mid, Y_q_mid);
% RSSI_q_Radar.mid = griddata(X, Y, RSSI_Radar, X_q_mid, Y_q_mid);

X_q_grid = linspace(1, 3, 5).*w;
Y_q_grid = linspace(1, -1, 5*N).*w;

[X_q_grid, Y_q_grid] = meshgrid(X_q_grid, Y_q_grid);

P_t_q.grid = griddata(X, Y, P_t, X_q_grid, Y_q_grid);
% RSSI_q.grid = griddata(X, Y, RSSI, X_q_grid, Y_q_grid);
% P_t_q_Radar.grid = griddata(X, Y, P_t_Radar, X_q_grid, Y_q_grid);
% RSSI_q_Radar.grid = griddata(X, Y, RSSI_Radar, X_q_grid, Y_q_grid);

%% Import Measurements.
% [G_R_T_NA, G_T_R_NA] = importS2P_2(TX.f, 'Data\Measurement');
% 
% P_a_t = -15;
% 
% P_t_NA.mid = G_R_T_NA.mid + P_a_t;
% P_t_NA.grid = G_R_T_NA.grid + P_a_t;
% RSSI_NA.mid = G_T_R_NA.mid + P_t_NA.mid + 10*log10(real(RX.Z)/Z_L);
% RSSI_NA.grid = G_T_R_NA.grid + P_t_NA.grid + 10*log10(real(RX.Z)/Z_L);
% 
% P_t_NA.mid = P_t_NA.mid(3:end);
% P_t_NA.grid = P_t_NA.grid(:, 2:end);
% RSSI_NA.mid = RSSI_NA.mid(3:end);
% RSSI_NA.grid = RSSI_NA.grid(:, 2:end);
% 
% X_NA_mid = linspace(1, 3, 9).*w;
% Y_NA_mid = 0;
% 
% X_NA_grid = linspace(1, 3, 5).*w;
% Y_NA_grid = linspace(1, -1, 5).*w;
% 
% [X_NA_grid, Y_NA_grid] = meshgrid(X_NA_grid, Y_NA_grid);

%% Import Simulations.
% [G_R_T_CST, G_T_R_CST] = importS2P_1(TX.f, 'Data\Simulation');
% 
% P_t_CST.mid = G_R_T_CST.mid + P_a_t ;
% P_t_CST.grid = G_R_T_CST.grid + P_a_t;
% RSSI_CST.mid = G_T_R_CST.mid + P_t_CST.mid + 10*log10(real(RX.Z)/Z_L);
% RSSI_CST.grid = G_T_R_CST.grid + P_t_CST.grid + 10*log10(real(RX.Z)/Z_L);
% 
% P_t_CST.mid = P_t_CST.mid(11:end);
% P_t_CST.grid = G_R_T_CST.grid + P_a_t;
% RSSI_CST.mid = RSSI_CST.mid(11:end);
% RSSI_CST.grid = G_T_R_CST.grid + P_t_CST.grid + 10*log10(real(RX.Z)/Z_L);
% 
% X_CST_mid = linspace(1, 3, 41).*w;
% Y_CST_mid = 0;
% 
% X_CST_grid = linspace(0.5, 3, 6).*w;
% Y_CST_grid = linspace(1, -1, 41).*w;
% 
% [X_CST_grid, Y_CST_grid] = meshgrid(X_CST_grid, Y_CST_grid);
% 
% X_CST_grid_new = linspace(1, 3, 5).*w;
% Y_CST_grid_new = linspace(1, -1, 21).*w;
% 
% [X_CST_grid_new, Y_CST_grid_new] = meshgrid(X_CST_grid_new, Y_CST_grid_new);
% 
% P_t_CST.grid = griddata(X_CST_grid, Y_CST_grid, P_t_CST.grid, X_CST_grid_new, Y_CST_grid_new);
% 
% X_CST_grid = X_CST_grid_new ;
% Y_CST_grid= Y_CST_grid_new;

%% Simulation.
figure('Position',  [500, 50, 700, 350])
hold on
plot(X_q_mid./w, P_t_q.mid, 'LineWidth', 3)  
% plot(X_q_mid./w, P_t_q_Radar.mid, ':', 'LineWidth', 3)
% plot(X_CST_mid./w, P_t_CST.mid, 'kx', 'LineWidth', 2, 'MarkerSize', 7.5)
hold off

% % legend('{\small Novel Method}', '{\small Friis Formula}', '{\small CST}', 'FontSize', 12, 'Interpreter', 'none')
% xlabel({'','{\small Distance from Reader [$\lambda$]}'}, 'Interpreter', 'none');
% xlim([1 3])
% xticks(1:.5:3)
% ylabel({'{\small Power on Tag [dBm]}',''}, 'Interpreter', 'none');
% ylim([-40 -20])
% yticks(-40 : 5 : -20)
% 
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
% gca_hold = gca;
% gca_hold.XRuler.TickLabelGapOffset = 15;
% gca_hold.YRuler.TickLabelGapOffset = 25;
% 
% grid on
% box on
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\Paper\SVG\Simulation\', filetag, '\0.svg');
% % saveas(gcf, fileaddress)
% 
for i = 1 : 5
    
    figure('Position',  [500, 50, 700, 350])
    hold on
    plot(Y_q_grid(:,1)./w, P_t_q.grid(:,i), 'LineWidth', 3)
%     plot(Y_q_grid(:,1)./w, P_t_q_Radar.grid(:,i), ':', 'LineWidth', 3)
%     plot(Y_CST_grid(:,i)./w, P_t_CST.grid(:,i), 'kx', 'LineWidth', 2, 'MarkerSize', 7.5)
    hold off

%     legend('{\small Novel Method}', '{\small Friis Formula}', '{\small CST Simulation}', 'Orientation', 'horizontal', 'Location', 'S', 'FontSize', 12, 'Interpreter', 'none')
    xlabel({'','{\small Misalignment with Reader [$\lambda$]}'}, 'Interpreter', 'none');
    xticks(-1:.5:1)
    ylabel({'{\small Power on Tag [dBm]}',''}, 'Interpreter', 'none');
    ylim([-75 -25])
%     yticks()
    
    set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
    gca_hold = gca;
    gca_hold.XRuler.TickLabelGapOffset = 15;
    gca_hold.YRuler.TickLabelGapOffset = 25;
    
    grid on
    box on
    
    filename = sprintf('%.1f', i/2);
    fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\Paper\SVG\Simulation\', filetag, '\', filename, '.svg');
%     saveas(gcf, fileaddress)
    
end

%% Measurement.
% figure('Position',  [500, 50, 700, 350])
% % yyaxis left
% hold on
% plot(X_q_mid./w, P_t_q.mid, 'LineWidth', 3)
% plot(X_NA_mid./w, P_t_NA.mid, 'ks', 'LineWidth', 2, 'MarkerSize', 7.5)
% hold off
% 
% xlabel({'','{\small Distance from Reader [$\lambda$]}'}, 'Interpreter', 'none');
% xlim([1 3])
% xticks(1:.5:3)
% ylabel({'{\small Magnitude [dBm]}',''}, 'Interpreter', 'none');
% ylim([-65 -25])
% yticks(-65 : 10 : -25)
% 
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
% gca_hold = gca;
% gca_hold.XRuler.TickLabelGapOffset = 15;
% gca_hold.YRuler.TickLabelGapOffset = 25;
% 
% grid on
% box on
% 
% % yyaxis right
% hold on
% plot(X_q_mid./w, RSSI_q.mid, ':', 'LineWidth', 3)
% plot(X_NA_mid./w, RSSI_NA.mid, 'kd', 'LineWidth', 2, 'MarkerSize', 7.5)
% hold off
% % 
% % ylabel({'','','{\small RSSI [dBm]}',''}, 'Interpreter', 'none');
% %     ylim([-65 -25])
% %     yticks(-65 : 10 : -25)
% % 
% % gca_hold = gca;
% % gca_hold.XRuler.TickLabelGapOffset = 15;
% % gca_hold.YAxis(2).TickLabelGapOffset = 15;
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\Paper\SVG\Measurement\', filetag, '\0.svg');
% % saveas(gcf, fileaddress)
% 
% for i = 1 : 5
%     
%     figure('Position',  [500, 50, 700, 350])
% %     yyaxis left
%     hold on
%     plot(Y_q_grid(:,1)./w, P_t_q.grid(:,i), 'LineWidth', 3)
%     plot(Y_NA_grid(:,1)./w, P_t_NA.grid(:,i), 'ks', 'LineWidth', 2, 'MarkerSize', 7.5)
%     hold off
% 
%     xlabel({'','{\small Misalignment with Reader [$\lambda$]}'}, 'Interpreter', 'none');
%     xticks(-1:.5:1)
%     ylabel({'{\small Magnitude [dBm]}',''}, 'Interpreter', 'none');
%     ylim([-75 -25])
%     yticks(-75 : 10 : -25)
%     
%     set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
%     gca_hold = gca;
%     gca_hold.XRuler.TickLabelGapOffset = 15;
%     gca_hold.YRuler.TickLabelGapOffset = 25;
%     
%     grid on
%     box on
%     
% %     yyaxis right
%     hold on
%     plot(Y_q_grid(:,1)./w, RSSI_q.grid(:,i), ':', 'LineWidth', 3)
%     plot(Y_NA_grid(:,1)./w, RSSI_NA.grid(:,i), 'kd', 'LineWidth', 2, 'MarkerSize', 7.5)
%     hold off
% 
% %     legend('{\small $P_{t}$ by Novel Method}', '{\small $P_{t}$ by Measurement}', '{\small RSSI by Novel Method}', '{\small RSSI by Measurement}', 'Orientation', 'horizontal', 'Location', 'S', 'FontSize', 12, 'Interpreter', 'none')
% %     ylabel({'','','{\small RSSI [dBm]}'}, 'Interpreter', 'none');
% %     ylim([-65 -25])
% %     yticks(-65 : 10 : -25)
% %     
% %     set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
% %     gca_hold = gca;
% %     gca_hold.XRuler.TickLabelGapOffset = 15;
% %     gca_hold.YAxis(2).TickLabelGapOffset = 15;
% %     
% %     grid on
% %     box on
%     
%     filename = sprintf('%.1f', i/2+0.5);
%     fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\Paper\SVG\Measurement\', filetag, '\', filename, '.svg');
% %     saveas(gcf, fileaddress)
%     
% end
% 
% %% IMS Simulation.
% figure('Position',  [500, 50, 700, 350])
% hold on
% plot(X_q_mid./w, P_t_q.mid, 'LineWidth', 3)  
% plot(X_q_mid./w, P_t_q_Radar.mid, 'LineWidth', 3)
% 
% idx = 1 : 41;
% X_CST_mid(rem(idx, 2) == 1);
% 
% plot(X_CST_mid(rem(idx, 2) == 1)./w, P_t_CST.mid(rem(idx, 2) == 1), 'kx', 'LineWidth', 3, 'MarkerSize', 10)
% hold off
% 
% % legend('{\small Novel Method}', '{\small Friis Formula}', '{\small CST}', 'FontSize', 12, 'Interpreter', 'none')
% title('y = 0')
% xlabel('x [\lambda]')
% xlim([1 3])
% xticks(1:.5:3)
% ylabel('Power [dB]')
% ylim([-40 -20])
% yticks(-40 : 5 : -20)
% 
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontName', 'Franklin Gothic Book', 'FontSize', 20)
% grid on
% box on
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\PPT\SVG\Simulation\', filetag, '\0.svg');
% saveas(gcf, fileaddress)
% 
% for i = 1 : 5
%     
%     figure('Position',  [500, 50, 700, 350])
%     hold on
%     plot(Y_q_grid(:,1)./w, P_t_q.grid(:,i), 'LineWidth', 3)
%     plot(Y_q_grid(:,1)./w, P_t_q_Radar.grid(:,i), 'LineWidth', 3)
%     plot(Y_CST_grid(:,i)./w, P_t_CST.grid(:,i), 'kx', 'LineWidth', 3, 'MarkerSize', 10)
%     hold off
% 
% %     legend('{\small Novel Method}', '{\small Friis Formula}', '{\small CST Simulation}', 'Orientation', 'horizontal', 'Location', 'S', 'FontSize', 12, 'Interpreter', 'none')
%     title(sprintf('x = %.0f\\lambda', i/2+0.5))
%     xlabel('y [\lambda]')
%     xticks(-1:.5:1)
%     ylabel('Power [dB]')
%     ylim([-45 -25])
% %     yticks()
%     
%     set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontName', 'Franklin Gothic Book', 'FontSize', 20)
%     
%     grid on
%     box on
%     
%     filename = sprintf('%.1f', i/2+0.5);
%     fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\PPT\SVG\Simulation\', filetag, '\', filename, '.svg');
%     saveas(gcf, fileaddress)
%     
% end

%% IMS Measurement.
% figure('Position',  [500, 50, 700, 350])
% % yyaxis left
% hold on
% plot(X_q_mid./w, P_t_q.mid, 'LineWidth', 3)
% plot(X_NA_mid./w, P_t_NA.mid, 'ks', 'LineWidth', 2.5, 'MarkerSize', 10)
% hold off
% 
% title('y = 0')
% xlabel('x [\lambda]')
% xlim([1 3])
% xticks(1:.5:3)
% ylabel('Power [dB]')
% ylim([-65 -25])
% yticks(-65 : 10 : -25)
% 
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontName', 'Franklin Gothic Book', 'FontSize', 20)
% grid on
% box on
% 
% % yyaxis right
% hold on
% plot(X_q_mid./w, RSSI_q.mid, 'LineWidth', 3)
% plot(X_NA_mid./w, RSSI_NA.mid, 'kd', 'LineWidth', 2.5, 'MarkerSize', 8.5)
% hold off
% % 
% % ylabel({'','','{\small RSSI [dBm]}',''}, 'Interpreter', 'none');
% %     ylim([-65 -25])
% %     yticks(-65 : 10 : -25)
% % 
% % gca_hold = gca;
% % gca_hold.XRuler.TickLabelGapOffset = 15;
% % gca_hold.YAxis(2).TickLabelGapOffset = 15;
% 
% % legend('1', '2', '3', '4')
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\PPT\SVG\Measurement\', filetag, '\0.svg');
% saveas(gcf, fileaddress)
% 
% for i = 1 : 5
%     
%     figure('Position',  [500, 50, 700, 350])
% %     yyaxis left
%     hold on
%     plot(Y_q_grid(:,1)./w, P_t_q.grid(:,i), 'LineWidth', 3)
%     plot(Y_NA_grid(:,1)./w, P_t_NA.grid(:,i), 'ks', 'LineWidth', 2.5, 'MarkerSize', 10)
%     hold off
%     
%     title(sprintf('x = %.0f\\lambda', i/2+0.5))
%     xlabel('y [\lambda]')
%     xticks(-1:.5:1)
%     ylabel('Power [dB]')
%     ylim([-75 -25])
%     yticks(-75 : 10 : -25)
%     
%     set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontName', 'Franklin Gothic Book', 'FontSize', 20)
%     
%     grid on
%     box on
%     
% %     yyaxis right
%     hold on
%     plot(Y_q_grid(:,1)./w, RSSI_q.grid(:,i), 'LineWidth', 3)
%     plot(Y_NA_grid(:,1)./w, RSSI_NA.grid(:,i), 'kd', 'LineWidth', 2.5, 'MarkerSize', 8.5)
%     hold off
% 
% %     legend('{\small $P_{t}$ by Novel Method}', '{\small $P_{t}$ by Measurement}', '{\small RSSI by Novel Method}', '{\small RSSI by Measurement}', 'Orientation', 'horizontal', 'Location', 'S', 'FontSize', 12, 'Interpreter', 'none')
% %     ylabel({'','','{\small RSSI [dBm]}'}, 'Interpreter', 'none');
% %     ylim([-65 -25])
% %     yticks(-65 : 10 : -25)
% %     
% %     set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
% %     gca_hold = gca;
% %     gca_hold.XRuler.TickLabelGapOffset = 15;
% %     gca_hold.YAxis(2).TickLabelGapOffset = 15;
% %     
% %     grid on
% %     box on
%     
%     filename = sprintf('%.1f', i/2+0.5);
%     fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\PPT\SVG\Measurement\', filetag, '\', filename, '.svg');
%     saveas(gcf, fileaddress)
%     
% end
% 
% %% Threshold, interpolate..
% P_t_c.mid = griddata(X, Y, P_t, X_NA_mid, Y_NA_mid);
% RSSI_c.mid = griddata(X, Y, RSSI, X_NA_mid, Y_NA_mid);
% P_t_Radar_c.mid = griddata(X, Y, P_t_Radar, X_NA_mid, Y_NA_mid);
% RSSI_Radar_c.mid = griddata(X, Y, RSSI_Radar, X_NA_mid, Y_NA_mid);
% 
% P_t_c.grid = griddata(X, Y, P_t, X_NA_grid, Y_NA_grid);
% RSSI_c.grid = griddata(X, Y, RSSI, X_NA_grid, Y_NA_grid);
% P_t_Radar_c.grid = griddata(X, Y, P_t_Radar, X_NA_grid, Y_NA_grid);
% RSSI_Radar_c.grid = griddata(X, Y, RSSI_Radar, X_NA_grid, Y_NA_grid);
% 
% %% Threshold.
% T_th = -34;
% 
% B_R_T.mid = P_t_c.mid >= T_th;
% B_R_T.grid = P_t_c.grid >= T_th;
% 
% B_R_T_Radar.mid = P_t_Radar_c.mid >= T_th;
% B_R_T_Radar.grid = P_t_Radar_c.grid >= T_th;
% 
% B_R_T_NA.mid = P_t_NA.mid >= T_th;
% B_R_T_NA.grid = P_t_NA.grid >= T_th;
% 
% Accuracy.Downlink.M_NA = (sum(~xor(B_R_T.mid, B_R_T_NA.mid), 'all') + sum(~xor(B_R_T.grid, B_R_T_NA.grid), 'all'))/41*100;
% Accuracy.Downlink.R_NA = (sum(~xor(B_R_T_Radar.mid, B_R_T_NA.mid), 'all') + sum(~xor(B_R_T_Radar.grid, B_R_T_NA.grid), 'all'))/41*100;
% 
% figure
% 
% hold on
% 
% R_polar = Rho_max;
% T_polar = linspace(0,2*pi);
% X_polar = R_polar.*cos(T_polar);
% Y_polar = R_polar.*sin(T_polar);
% 
% plot(X_polar, Y_polar, 'k', 'LineWidth', 1)
% 
% for i = 1 : 2
%     
%     plot(i.*X_polar./3, i.*Y_polar./3, 'Color', [.5 .5 .5], 'LineStyle', ':')
% 
% end
% 
% R_polar = Rho_max.*linspace(-1, 1);
% T_polar = 0:30:360;
% X_polar = R_polar.*cos(deg2rad(T_polar))';
% Y_polar = R_polar.*sin(deg2rad(T_polar))';
% for i = 1 : length(T_polar)
%    
%     plot(X_polar(i,:), Y_polar(i,:), 'Color', [.5 .5 .5], 'LineStyle', ':')
%     
% end
% 
% contour(X, Y, P_t, linspace(-70, 30, 11) + T_th, 'LineWidth', .5, 'LineColor', 'k', 'ShowText', 'off');
% [~, h2] = contourf(X, Y, P_t, [T_th T_th], 'LineStyle', ':','LineWidth', 2, 'Color', 'k', 'ShowText', 'off');
% colormap(winter)
% h3 = plot(X_NA_grid(~B_R_T_NA.grid), Y_NA_grid(~B_R_T_NA.grid), 'x', 'LineWidth', .5,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k');
% h4 = plot(X_NA_grid(B_R_T_NA.grid), Y_NA_grid(B_R_T_NA.grid), 'o', 'LineWidth', 1.5,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k');
% plot(X_NA_mid(~B_R_T_NA.mid), zeros(size(X_NA_mid(~B_R_T_NA.mid))), 'x', 'LineWidth', .5,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k')
% plot(X_NA_mid(B_R_T_NA.mid), zeros(size(X_NA_mid(B_R_T_NA.mid))), 'o', 'LineWidth', 1.5,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k')
% 
% axis equal
% caxis([10*T_th T_th+5]) 
% set(gca, 'visible', 'off')
% 
% % legend([h2 h3 h4], '{\small Simulated Tag Range}', '{\small Measured Tag Off}', '{\small Measured Tag On}', 'Position', [.775 .125 .2 0], 'FontSize', 12, 'Interpreter', 'none')
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\SVG\Threshold\', filetag, '\matlab.svg');
% % saveas(gcf, fileaddress)
% 
% figure
% 
% hold on
% 
% R_polar = Rho_max;
% T_polar = linspace(0,2*pi);
% X_polar = R_polar.*cos(T_polar);
% Y_polar = R_polar.*sin(T_polar);
% 
% plot(X_polar, Y_polar, 'k', 'LineWidth', 1)
% 
% for i = 1 : 2
%     
%     plot(i.*X_polar./3, i.*Y_polar./3, 'Color', [.5 .5 .5], 'LineStyle', ':')
% 
% end
% 
% R_polar = Rho_max.*linspace(-1, 1);
% T_polar = 0:30:360;
% X_polar = R_polar.*cos(deg2rad(T_polar))';
% Y_polar = R_polar.*sin(deg2rad(T_polar))';
% for i = 1 : length(T_polar)
%    
%     plot(X_polar(i,:), Y_polar(i,:), 'Color', [.5 .5 .5], 'LineStyle', ':')
%     
% end
% 
% contour(X, Y, P_t_Radar, linspace(-70, 30, 11) + T_th, 'LineWidth', .5, 'LineColor', 'k', 'ShowText', 'off');
% [~, h2] = contourf(X, Y, P_t_Radar, [T_th T_th], 'LineStyle', ':','LineWidth', 2, 'Color', 'k', 'ShowText', 'off');
% colormap(jet)
% h3 = plot(X_NA_grid(~B_R_T_NA.grid), Y_NA_grid(~B_R_T_NA.grid), 'x', 'LineWidth', .5,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k');
% h4 = plot(X_NA_grid(B_R_T_NA.grid), Y_NA_grid(B_R_T_NA.grid), 'o', 'LineWidth', 1.5,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k');
% plot(X_NA_mid(~B_R_T_NA.mid), zeros(size(X_NA_mid(~B_R_T_NA.mid))), 'x', 'LineWidth', .5,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k')
% plot(X_NA_mid(B_R_T_NA.mid), zeros(size(X_NA_mid(B_R_T_NA.mid))), 'o', 'LineWidth', 1.5,...
%     'MarkerSize',3,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k')
% 
% axis equal
% caxis([1.3*T_th T_th+5]) 
% set(gca, 'visible', 'off')
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\SVG\Threshold\', filetag, '\radar.svg');
% % saveas(gcf, fileaddress)
% 
% %% IMS.
% T_th = -34;
% 
% figure('Position',  [500, 50, 770, 385])
% 
% hold on
% contour(X./w, Y./w, P_t, [T_th T_th], 'LineStyle', ':','LineWidth', 3, 'Color', 'k', 'ShowText', 'off');
% contour(X./w, Y./w, P_t_Radar, [T_th T_th], 'LineStyle', ':','LineWidth', 3, 'Color', 'k', 'ShowText', 'off');
% 
% h3 = plot(X_NA_grid(~B_R_T_NA.grid)./w, Y_NA_grid(~B_R_T_NA.grid)./w, 'x', 'LineWidth', .5,...
%     'MarkerSize',7.5,...
%     'MarkerEdgeColor','k');
% h4 = plot(X_NA_grid(B_R_T_NA.grid)./w, Y_NA_grid(B_R_T_NA.grid)./w, 'o', 'LineWidth', 1.5,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k');
% plot(X_NA_mid(~B_R_T_NA.mid)./w, zeros(size(X_NA_mid(~B_R_T_NA.mid))), 'x', 'LineWidth', .5,...
%     'MarkerSize',7.5,...
%     'MarkerEdgeColor','k')
% plot(X_NA_mid(B_R_T_NA.mid)./w, zeros(size(X_NA_mid(B_R_T_NA.mid))), 'o', 'LineWidth', 1.5,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','k')
% 
% corr = zeros(size(TX.Layout.Points));
% corr(:,1) = .5;
% 
% % Layout = triangulation(TX.Layout.ConnectivityList, TX.Layout.Points./1000-corr);
% % trisurf(Layout)
% 
% axis equal
% xlim([-1 3.5])
% ylim([-1.5 1.5])
% xlabel('x [\lambda]');
% ylabel('y [\lambda]');
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'FontName', 'Franklin Gothic Book', 'FontSize', 20)
% 
% grid on
% box on
% 
% % legend([h2 h3 h4], '{\small Simulated Tag Range}', '{\small Measured Tag Off}', '{\small Measured Tag On}', 'Position', [.775 .125 .2 0], 'FontSize', 12, 'Interpreter', 'none')
% 
% fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\IMS\PPT\SVG\Threshold\', filetag, '\threshold.svg');
% saveas(gcf, fileaddress)
% 
% %%
% close all