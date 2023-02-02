clc
clear all
close all

f = 866E6;
w = physconst('Lightspeed')/f;
k = 2*pi/w;
r = linspace(0, 15/k, 10000);

x = k.*r;

for l = 0 : 2
    
    j = -besselj(l-1, x);
    y = bessely(l, x);
    e(l+1, :) = j - y;
    
    figure('Position',  [500, 50, 700, 350])
    hold on
    plot(x, j, 'LineWidth', 3)
    plot(x, y, ':', 'LineWidth', 3)
    plot(ones(1,10).*l, linspace(-3, 1, 10), 'k:', 'LineWidth', 1)
    hold off
    
    legend('{\small Modified}', '{\small Original}', 'Location', 'SE', 'FontSize', 12, 'Interpreter', 'none')
    xlabel({'','{\small Argument $kr$}'}, 'Interpreter', 'none')
    ylabel({'{\small Neumann Function}',''}, 'Interpreter', 'none')
    ylim([-1.5 1])
    
    set(gca, 'XTick', unique([xticks, l]), 'GridLineStyle', ':', 'GridAlpha', 0.5)
    gca_hold = gca;
    gca_hold.XRuler.TickLabelGapOffset = 15;
    gca_hold.YRuler.TickLabelGapOffset = 25;
    
    grid on
    box on
    
    filename = sprintf('bessel_%d', l);
    fileaddress = strcat('C:\Users\Administrator\OneDrive - UGent\RFID\Paper\SVG\Bessel\', filename, '.svg');
%     saveas(gcf, fileaddress)
    
end