function [lon, lat] = groundTrack (state, theta_G0, t, w_E, t_0)
%
% PROTOTYPE:
%  [alpha, delta, lon, lat] = groundTrack (state, theta_G0, t, w_E, t_0);
%
% DESCRIPTION:
% Function to obtain groundtrack variables from cartesian coordinates
%
% INPUT:
%   state [?,6]     position and velocity matrix [km, km/s]
%   theta_G0 [1]    initial Greenwich Meridian Longitude [rad]
%   t [?,1]         time vector from ODE [s]
%   w_E [1]         Earth angular velocity [rad/s]
%   t_0 [1]         initial value of time [s]
%
% OUTPUT:
%   lon [?,1]       longitude [deg]
%   lat [?,1]       latitude [deg]
%
% FUNCTIONS CALLED:
% (none)
%
% AUTHORS: 
%   Giuseppe Brentino, Virginia di Biagio Missaglia, Nicolò
%   Galletta, Roberto Pistone Nascone, 2022
% -------------------------------------------------------------------------

r = sqrt(state(:, 1).^2 + state(:, 2).^2 + state(:, 3).^2);

delta = asin(state(:, 3)./r(:));            % declination [rad]
alpha = atan2(state(:, 2), state(:, 1));    % right ascension [rad]

theta_G = theta_G0 + w_E.*(t(:) - t_0);
lon = alpha(:) - theta_G(:);
lat = delta(:);

lon(:)= wrapToPi(lon(:));

lat = rad2deg(lat);
lon = rad2deg(lon);
