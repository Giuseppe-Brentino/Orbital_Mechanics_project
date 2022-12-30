function animated_plot(X,Y,Z, n_periods, n_point, inclination)
%
% PROTOTYPE:
% animated_plot(X,Y,Z, n_periods, n_point, inclination)
%
% DESCRIPTION:
% The function creates an animated plot of the orbit evolutions over an
% integer multiple of time periods. This function works fine only when
% orbit propagation is performed with fixed time step for each period.
% Variations over time of orbital periods are assumed to be small.
%
% INPUT:
%   X [?,1]         x coordinates of S/C position vector in cartesian rf [km]
%   X [?,1]         y coordinates of S/C position vector in cartesian rf [km]
%   X [?,1]         z coordinates of S/C position vector in cartesian rf [km]
%   n_periods [1]   number of periods of the plot [-]
%   n_point [1]     number of time interval for each period [-]
%   inclination [1] initial inclination angle [deg] (default is 0)
%
%   Note: n_periods and n_point must be positive integer numbers,
%   inclination can be omitted
%
% FUNCTIONS CALLED:
% plot_terra.m
%
% AUTHORS: 
%   Giuseppe Brentino, Virginia di Biagio Missaglia, Nicolò
%   Galletta, Roberto Pistone Nascone, 2022
% -------------------------------------------------------------------------


suplim = max([max(X), max(Y), max(Z)])+2000;
inflim = min([min(X), min(Y), min(Z)])-2000;

if nargin == 5
    inclination = 0;
end

figure()
plot_terra
set(gca, 'color', '#010020')
axis equal
hold on

h = animatedline('Color', '#00d400', 'LineWidth', 2);
axis([inflim suplim inflim suplim inflim suplim])
xlabel('r_x [km]')
ylabel('r_y [km]')
zlabel('r_z [km]')

for i=1:n_periods

    clearpoints(h)

    x = [ X( (n_point*i-n_point + 1): n_point*i)];
    y = [ Y( (n_point*i-n_point + 1): n_point*i)];
    z = [ Z( (n_point*i-n_point + 1): n_point*i)];

    addpoints(h, x, y, z)
    view(180, 90+inclination)
    title(['orbit evolution at time ' num2str(i) 'T'], 'Color', '#fdf6e4', ...
        'FontSize', 8)
    set(get(gca,'title'),'Position',[1/2*suplim inflim+1000 1/2*suplim])

    movievector(i) = getframe;

end

% save the animated plot
txt = input('Do you want to save the animated plot? yes/no \n', 's');

if txt == 'yes'
    
    myWriter = VideoWriter('orbit_animatedplot');
    open(myWriter);
    writeVideo(myWriter, movievector);
    close(myWriter)

end