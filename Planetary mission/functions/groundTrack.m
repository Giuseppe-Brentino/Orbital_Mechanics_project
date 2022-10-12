function [alpha, delta, lon, lat] = groundTrack_2 (state, theta_G0, t, w_E, settings, t_0)
%
% PROTOTYPE
% [alpha, delta, longitude, latitude] = graoundTrack(s_0, theta_G0, tspan, w_E, settings, t_0)
% Function to obtain groundtrack variables from chartesian coordinates
%
% ------------------------------INPUT-------------------------------------
% state = position and velocity matrix [?x6]
% theta_G0 = Greenwich Meridian Longitude [1x1]
% t = time vector from ODE [?x1] [s]
% w_E = Earth angular velocity [1x1] [rad/s]
% settings = settings to evaluate mu
% t_0 = initial value of time [1x1] [s]
%
% ------------------------------OUTPUT------------------------------------
% alpha =  right ascension [1x1] [deg]
% delta =  declination [1x1] [deg]
% longitude =  [1x1] [deg]
% latitude = [1x1] [deg]


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




















