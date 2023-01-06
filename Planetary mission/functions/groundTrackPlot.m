function groundTrackPlot(figureName)

% DESCRIPTION:
% Call this function in order to plot Earth as background of the
% groundtrack plots
%
% USE:
% Just call the function
%
% EARTH's IMAGE FROM: 
% Reto Stöckli, NASA Earth Observatory, Blue Marble Next Generation
% website: https://visibleearth.nasa.gov
%--------------------------------------------------------------------------

if nargin == 0
    figure()
else
    figure(figureName)
end

hold on;
grid on;
photo = imread('Earth.png');
image('XData', [-180, 180], 'YData', [-90, 90], 'CData', flipud(photo));
