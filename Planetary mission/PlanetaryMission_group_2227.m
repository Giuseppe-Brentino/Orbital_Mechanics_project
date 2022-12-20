clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');


% data
a = 1.9539*1e4;             % semi-major axis [km]
e = 0.6607;                 % eccentricity [-]
i = 10.4278;                % inclination [deg]
OM = 60;                    % RAAN [deg]
om = 30;                    % argument of pericentre [deg]
theta = 0;                  % true anomaly [deg]

GT_ratio = 13/4;            % satellite's revs wrt Earth's revs [-]
c_d = 2.1;                  % drag coefficient [-]
Area_mass = 0.0042;         % ratio between reference area and mass [m^2/kg]

settings.mu = astroConstants(13);               % Earth's gravitational constant [km^3/s^2]
settings.J2E = astroConstants(9);               % Earth oblateness parameter [-]
settings.RE = astroConstants(23);               % Earth radius [km]


%% nominal GT 

T = 2*pi*sqrt(a^3/settings.mu);             % Period of 1 orbit [s]
T_24 = 24*3600;                             % 1 day [s] 
T_10 = 10*T_24;                             % 10 days [s]
%T/T_24;

w_E = deg2rad(15.04)/3600;                                % Earth's angular velocity [rad/s]
[r0, v0] = kep2car(a, e, i, OM, om, theta, settings.mu);
s0 = [r0; v0];                              % initial state vector

options = odeset('RelTol', 1e-10,'AbsTol',1e-11 );
[t, Y] = ode113(@pert_tbp,[0 T], s0, options, settings);

t0 = 0;
theta_G0 = 0;

[alpha, delta, lon, lat] = groundTrack (Y, theta_G0, t, w_E, t0);

for j = 2:length(lon)
    if lon(j-1)*lon(j) < 0 && abs(lon(j)) > 10
        lon = [lon(1:j-1); NaN; lon(j:end)];
        lat = [lat(1:j-1); NaN; lat(j:end)];
   end

    if lat(j-1)*lat(j) < 0 && abs(lat(j)) > 10
        lon = [lon(1:j-1); NaN; lon(j:end)];
        lat = [lat(1:j-1); NaN; lat(j:end)];
    end

end

% TOGLIERE I MARGINI BIANCHI DAL PLOT DELLA GT
groundTrackPlot
% xlim([-180 180]);
% ylim([-90 90]);
plot(lon, lat, 'r','LineWidth',2);
plot(lon(1), lat(1), 'o', 'LineWidth', 5);
plot(lon(end), lat(end), '*', 'LineWidth', 5);


%% Modified GT

T_mod = (2*pi/w_E)/GT_ratio;
a_mod = ( settings.mu*(T_mod/(2*pi))^2 )^(1/3);

[r_mod, v_mod] = kep2car(a_mod, e, i, OM, om, theta, settings.mu);
s_mod = [r_mod; v_mod];                              % initial modified state vector

[t1, Y1] = ode113(@pert_tbp,[0 T_mod], s_mod, options, settings);

[alpha_mod, delta_mod, lon_mod, lat_mod] = groundTrack (Y1, theta_G0, t1, w_E, t0);

for j = 2:length(lon_mod)
    if lon_mod(j-1)*lon_mod(j) < 0 && abs(lon_mod(j)) > 10
        lon_mod = [lon_mod(1:j-1); NaN; lon_mod(j:end)];
        lat_mod = [lat_mod(1:j-1); NaN; lat_mod(j:end)];
   end

    if lat_mod(j-1)*lat_mod(j) < 0 && abs(lat_mod(j)) > 10
        lon_mod = [lon_mod(1:j-1); NaN; lon_mod(j:end)];
        lat_mod = [lat_mod(1:j-1); NaN; lat_mod(j:end)];
    end

end

% TOGLIERE I MARGINI BIANCHI DAL PLOT DELLA GT
groundTrackPlot
% xlim([-180 180]);
% ylim([-90 90]);
plot(lon_mod, lat_mod, 'r','LineWidth',2);
plot(lon, lat, 'g','LineWidth',2);
plot(lon(1), lat(1), 'o', 'LineWidth', 5);
plot(lon(end), lat(end), '*', 'LineWidth', 5);
plot(lon_mod(1), lat_mod(1), 'o', 'LineWidth', 5);
plot(lon_mod(end), lat_mod(end), '*', 'LineWidth', 5);


%% Perturbated GT

% r_norm = [];
% for k = 1:size(Y, 1)
%     r_norm(k) = norm(Y(k, 1:3));
% end
% 
% h = r_norm' - settings.RE;







