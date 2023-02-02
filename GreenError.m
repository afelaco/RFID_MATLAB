clc
clear all
close all

f = 866E6;
w = physconst('Lightspeed')/f;
k = 2*pi/w;
r = linspace(0, 3, 1000).*w;

for i = 1 : length(r)
    
    I_1 = (1 - 1i/(k*r(i)) - 1/(k*r(i))^2).*eye(3);
    I_2 = zeros(3,3);
    I_2(1,1) = 1;
    I_2 = (1 - 3i/(k*r(i)) - 3/(k*r(i))^2).*I_2;

    D = I_1 - I_2;
   
    A(i) = abs(D(2,2));
    b(i) = rad2deg(angle(D(2,2)));
    
end

%%
figure('Position',  [500, 50, 700, 350])
plot(r./w, abs((1-A)).*100, 'LineWidth', 3)  

% legend('{\small Novel Method}', '{\small Friis Formula}', '{\small CST}', 'FontSize', 12, 'Interpreter', 'none')
xlabel({'','{\small Distance [$\lambda$]}'}, 'Interpreter', 'none');
% xlim([0 5])
xticks(0:.5:3)
ylabel({'{\small Amplitude Error [%]}',''}, 'Interpreter', 'none');
ylim([0 20])
% yticks(-40 : 5 : -20)

set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5)
gca_hold = gca;
gca_hold.XRuler.TickLabelGapOffset = 15;
gca_hold.YRuler.TickLabelGapOffset = 25;

grid on
box on

fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\RFID\Paper\SVG\Green\amplitude.svg');
% saveas(gcf, fileaddress)

%%
figure('Position',  [500, 50, 700, 350])
plot(r./w, -b, 'LineWidth', 3, 'Color', [0.8500, 0.3250, 0.0980])  

% legend('{\small Novel Method}', '{\small Friis Formula}', '{\small CST}', 'FontSize', 12, 'Interpreter', 'none')
xlabel({'','{\small Distance [$\lambda$]}'}, 'Interpreter', 'none');
% xlim([0 5])
xticks(0:.5:3)
ylabel({'{\small Phase Error}',''}, 'Interpreter', 'none');
% ylim([-20 100])
yticks(0 : 45 : 180)
yticklabels({'0','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','{\small $\pi$}'})

set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.5, 'TickLabelInterpreter','none')
gca_hold = gca;
gca_hold.XRuler.TickLabelGapOffset = 15;
gca_hold.YRuler.TickLabelGapOffset = 25;

grid on
box on

fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\RFID\Paper\SVG\Green\phase.svg');
saveas(gcf, fileaddress)