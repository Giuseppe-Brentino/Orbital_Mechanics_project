function movievector = animated_plot(X,Y,Z, T, n_periods, p_span, n_point, pov, Eq_Plane)
%
% PROTOTYPE:
% movievector = animated_plot(X,Y,Z, T, n_periods, p_span, n_point, pov, Eq_Plane)
%
% DESCRIPTION:
% The function creates an animated plot of the orbit evolutions over an
% integer multiple of time periods. This function works fine only when
% orbit propagation is performed with fixed time step for each period.
% Variations over time of orbital periods are assumed to be small.
%
% INPUT:
%   X [?,1]         x coordinates of S/C position vector in cartesian rf [km]
%   Y [?,1]         y coordinates of S/C position vector in cartesian rf [km]
%   Z [?,1]         z coordinates of S/C position vector in cartesian rf [km]
%   T [1]           orbit initial period [s]
%   n_periods [1]   number of periods of the plot [-]
%   p_span [1]      periods for each frame [-]
%   n_point [1]     number of time interval for each period [-]
%   pov [2]         point of view: azimuth and elevation [deg] 
%   Eq_Plane        plot equatorial plane (true or false)
%
%   Note: n_periods, p_span and n_point must be positive integer numbers,
%   pov and Eq_Plane can be omitted (default: view=[0,0], Eq_Plane=false)
%
% OUTPUT:
%   movievector     struct containing all the frames of the animated plot
%
% FUNCTIONS CALLED:
% plot_terra.m
%
% AUTHORS: 
%   Giuseppe Brentino, Virginia di Biagio Missaglia, Nicolò
%   Galletta, Roberto Pistone Nascone, 2022
% -------------------------------------------------------------------------

suplim = max([max(X), max(Y), max(Z)]) +1000;
inflim = min([min(X), min(Y), min(Z)]) -1000;

k=max(suplim, abs(inflim));

if nargin == 7
    pov = [0,0];
    Eq_Plane=false;
end

az=pov(1);
el=pov(2);

figure('Units','normalized','OuterPosition', [0 0 1 1])
plot_terra
set(gca, 'color', '#010020')
axis equal

if Eq_Plane
patch([1 -1 -1 1]*k, [1 1 -1 -1]*k, [0 0 0 0], 'b', 'FaceAlpha', 0.8)
end

h = animatedline('Color', '#00d400', 'LineWidth', 2);
axis([inflim suplim inflim suplim -7000 7000])

k=1;

for i=1:p_span:n_periods

    clearpoints(h)

    x = [ X( (n_point*i-n_point + 1): n_point*i)];
    y = [ Y( (n_point*i-n_point + 1): n_point*i)];
    z = [ Z( (n_point*i-n_point + 1): n_point*i)];

    year = i*T/(365.25*24*3600);

    addpoints(h, x, y, z)
    view(az, el)
    title(['elapsed years: ' num2str(year)], ...
        'Color', '#fdf6e4', 'FontSize', 11)
    subtitle(['elapsed periods: ' num2str(i)], ...
        'Color', '#fdf6e4', 'FontSize', 11)
    set(get(gca,'title'),'Position', [1/2*suplim 1/2*suplim 1/2*suplim])
    set(get(gca,'subtitle'),'Position', [1/3*suplim 1/3*suplim+5000 1/3*suplim+3000])
    set(gca,'XColor', 'none','YColor','none', 'ZColor', 'none')

    movievector(k) = getframe;
    k=k+1;

end

% save the animated plot
txt = input('Do you want to save the animated plot? yes/no \n', 's');

if strcmp(txt,'yes')
    
    myWriter = VideoWriter('orbitAnimata');
    myWriter.FrameRate = (k-1)/20; % (20 seconds is the movie duration)

    open(myWriter);
    writeVideo(myWriter, movievector);
    close(myWriter)

end