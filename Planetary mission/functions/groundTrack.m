function [alpha, delta, lon, lat] = groundTrack (state, theta_G0, t, w_E, t_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function to obtain groundtrack variables from cartesian coordinates
%
% INPUT:
%   state = position and velocity matrix [?x6]
%   theta_G0 = Initial Greenwich Meridian Longitude [1x1]
%   t = time vector from ODE [?x1] [s]
%   w_E = Earth angular velocity [1x1] [rad/s]
%   t_0 = initial value of time [1x1] [s]
%
% OUTPUT:
%   alpha =  right ascension [?x1] [deg]
%   delta =  declination [?x1] [deg]
%   longitude =  [?x1] [deg]
%   latitude = [?x1] [deg]
%
% Contributors: 
%   Giuseppe Brentino, Virginia di Biagio Missaglia, Nicolò
%   Galletta, Roberto Pistone Nascone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta = zeros(length(t), 1);
alpha = zeros(length(t), 1);
lon = zeros(length(t), 1);
lat = zeros(length(t), 1);
theta_G = zeros(length(t), 1);

r = sqrt(state(:, 1).^2 + state(:, 2).^2 + state(:, 3).^2);

delta = asin(state(:, 3)./r(:));
alpha = atan2(state(:, 2), state(:, 1));

theta_G = theta_G0 + w_E.*(t(:) - t_0);
lon = alpha(:) - theta_G(:);
lat = delta(:);

lon(:)= wrapToPi(lon(:));

lat = rad2deg(lat);
lon = rad2deg(lon);




















