figure()
hold on;
grid on;
a = imread('Earth.png');
image('XData', [-180, 180], 'YData', [-90, 90], 'CData', flipud(a));
