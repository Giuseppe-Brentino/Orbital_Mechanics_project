clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');


% data
a = 1.9539*1e4;             % semi-major axis [km]
e = 0.6607;                 % eccentricity [-]
i = deg2rad(10.4278);       % inclination [rad]
OM = deg2rad(60);           % RAAN [rad]
om = deg2rad(30);           % argument of pericentre [rad]
theta = 0;                  % true anomaly [rad]

GT_ratio = 13/4;            % satellite's revs wrt Earth's revs [-]
drag.c_d = 2.1;             % drag coefficient [-]
drag.Area_mass = 0.0042;    % ratio between reference area and mass [m^2/kg]

settings.mu = astroConstants(13);           % Earth's gravitational constant [km^3/s^2]
settings.J2E = astroConstants(9);           % Earth oblateness parameter [-]
settings.RE = astroConstants(23);           % Earth radius [km]
settings.w_E = deg2rad(15.04)/3600;         % Earth's angular velocity [rad/s]

T = 2*pi*sqrt(a^3/settings.mu);             % Period of 1 orbit [s]
T_24 = 24*3600;                             % 1 day [s] 
T_10 = 10*T_24;                             % 10 days [s]
T_year = 365.25*T_24;                       % 1 year


%% nominal GT 

t_span = 300;
t_sample = 0:T/t_span:T;

%T/T_24;
[r0, v0] = kep2car(a, e, i, OM, om, theta, settings.mu);
s0 = [r0; v0];                              % initial state vector
settings.perturbations = false;
options = odeset('RelTol', 1e-10,'AbsTol',1e-11 );
[t, Y] = ode113(@pert_tbp,t_sample, s0, options, settings, drag);

t0 = 0;
theta_G0 = 0;
[alpha, delta, lon, lat] = groundTrack (Y, theta_G0, t, settings.w_E, t0);

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

groundTrackPlot
title('Nominal ground track')
xlim([-180 180]);
ylim([-90 90]);
xlabel('longitude [deg]')
ylabel('latitude [deg]')
plot(lon, lat, 'r','LineWidth',1);
plot(lon(1), lat(1), 'o', 'LineWidth', 1);
plot(lon(end), lat(end), '*', 'LineWidth', 1);
legend('groundtrack', 'initial position', 'final position')


%% Modified GT
T_mod = (2*pi/settings.w_E)/GT_ratio;
a_mod = ( settings.mu*(T_mod/(2*pi))^2 )^(1/3);

t_sample_mod = 0:T_mod/t_span:T_mod;

[r_mod, v_mod] = kep2car(a_mod, e, i, OM, om, theta, settings.mu);
s_mod = [r_mod; v_mod];                              % initial modified state vector

settings.perturbations = false;
[t_mod, Y_mod] = ode113(@pert_tbp, t_sample_mod, s_mod, options, settings, drag);

[alpha_mod, delta_mod, lon_mod, lat_mod] = groundTrack (Y_mod, theta_G0, t_mod, settings.w_E, t0);

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

groundTrackPlot
title('Modified ground track')
xlim([-180 180]);
ylim([-90 90]);
xlabel('longitude [deg]')
ylabel('latitude [deg]')
plot(lon_mod, lat_mod, 'r','LineWidth',1);
plot(lon_mod(1), lat_mod(1), 'or', 'LineWidth', 1);
plot(lon_mod(end), lat_mod(end), '*r', 'LineWidth', 1);
legend('modified ground track', 'modified orbit initial position', ...
    'modified orbit final position')

groundTrackPlot
title('Nominal and Modified ground track')
xlim([-180 180]);
ylim([-90 90]);
xlabel('longitude [deg]')
ylabel('latitude [deg]')
plot(lon_mod, lat_mod, 'r','LineWidth',1);
plot(lon_mod(1), lat_mod(1), 'or', 'LineWidth', 1);
plot(lon_mod(end), lat_mod(end), '*r', 'LineWidth', 1);
plot(lon, lat, 'g','LineWidth',1);
plot(lon(1), lat(1), 'og', 'LineWidth', 1);
plot(lon(end), lat(end), '*g', 'LineWidth', 1);
legend('modified ground track', 'modified orbit initial position', ...
    'modified orbit final position', 'nominal ground track', ...
    'nominal orbit initial position', 'nominal orbit final position')


%% Perturbed GT

settings.perturbations = true;
[t_pert, Y_pert] = ode113(@pert_tbp,t_sample, s0, options, settings, drag);

figure()
title('Orbit plots')
hold on;
grid on;
axis equal;
xlabel('r_x [km]')
ylabel('r_y [km]')
zlabel('r_z [km]')
plot_terra;
plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'DisplayName', 'Nominal orbit');
plot3(Y_pert(:, 1), Y_pert(:, 2), Y_pert(:, 3), 'DisplayName', 'Perturbed orbit');

[alpha_p, delta_p, lon_p, lat_p] = groundTrack (Y_pert, theta_G0, t_pert, settings.w_E, t0);

for j = 2:length(lon_p)
    if lon_p(j-1)*lon_p(j) < 0 && abs(lon_p(j)) > 10
        lon_p = [lon_p(1:j-1); NaN; lon_p(j:end)];
        lat_p = [lat_p(1:j-1); NaN; lat_p(j:end)];
   end

    if lat_p(j-1)*lat_p(j) < 0 && abs(lat_p(j)) > 10
        lon_p = [lon_p(1:j-1); NaN; lon_p(j:end)];
        lat_p = [lat_p(1:j-1); NaN; lat_p(j:end)];
    end

end


groundTrackPlot
title('Perturbed ground track')
xlim([-180 180]);
ylim([-90 90]);
xlabel('longitude [deg]')
ylabel('latitude [deg]')
plot(lon_p, lat_p, 'y','LineWidth',1);
plot(lon_p(1), lat_p(1), 'o', 'LineWidth', 1);
plot(lon_p(end), lat_p(end), '*', 'LineWidth', 1);
legend('perturbed gt', 'initial position', 'final position')

%% Perturbed and Modified GT

settings.perturbations = true;   
[t_p_mod, Y_p_mod] = ode113(@pert_tbp,t_sample_mod, s_mod, options, settings, drag);

[alpha_p_mod, delta_p_mod, lon_p_mod, lat_p_mod] = groundTrack (Y_p_mod, ...
    theta_G0, t_p_mod, settings.w_E, t0);

for j = 2:length(lon_p_mod)
    if lon_p_mod(j-1)*lon_p_mod(j) < 0 && abs(lon_p_mod(j)) > 10
        lon_p_mod = [lon_p_mod(1:j-1); NaN; lon_p_mod(j:end)];
        lat_p_mod = [lat_p_mod(1:j-1); NaN; lat_p_mod(j:end)];
   end

    if lat_p_mod(j-1)*lat_p_mod(j) < 0 && abs(lat_p_mod(j)) > 10
        lon_p_mod = [lon_p_mod(1:j-1); NaN; lon_p_mod(j:end)];
        lat_p_mod = [lat_p_mod(1:j-1); NaN; lat_p_mod(j:end)];
    end

end

groundTrackPlot
title('Perturbed and modified ground track')
xlim([-180 180]);
ylim([-90 90]);
xlabel('longitude [deg]')
ylabel('latitude [deg]')
plot(lon_p_mod, lat_p_mod, 'y','LineWidth',1);
plot(lon_p_mod(1), lat_p_mod(1), 'o', 'LineWidth', 1);
plot(lon_p_mod(end), lat_p_mod(end), '*', 'LineWidth', 1);
legend('perturbed and modified gt', 'initial position', 'final position')


%% Orbit propagation with Gauss's Eqns

kep0 = [a e i OM om theta];
[t_G, kep_G] = ode113(@kep_pert, t_sample, kep0, options, settings, drag);

Y_G = zeros(size(kep_G,1), 6);
for j = 1:size(kep_G, 1)
    [Y_G(j, 1:3), Y_G(j, 4:6)] = kep2car(kep_G(j, 1), kep_G(j, 2), kep_G(j, 3), kep_G(j, 4), kep_G(j, 5), kep_G(j, 6));

end

figure()
title('Orbit plots')
hold on;
grid on;
axis equal;
xlabel('r_x [km]')
ylabel('r_y [km]')
zlabel('r_z [km]')
plot_terra;
plot3(Y_pert(:, 1), Y_pert(:, 2), Y_pert(:, 3));
plot3(Y_G(:, 1), Y_G(:, 2), Y_G(:, 3));
legend('Cartesian propagation', 'Gauss propagation');

%% History of Keplerian elements

kep_pert = zeros(size(Y_pert,1),6);
for j = 1:size(Y_pert,1)

    r = Y_pert(j, 1:3);
    v = Y_pert(j, 4:6);

    [a, e, i, OM, om, th] = car2kep(r, v, settings.mu);
    kep_pert(j, :) = [a, e, rad2deg(i), rad2deg(OM), rad2deg(om), th];

end

kep_pert(:, 6) = rad2deg(unwrap(kep_pert(:, 6)));

for j = 1:length(kep_G)

    kep_G(j, 3:5) = [rad2deg(kep_G(j,3)), rad2deg(kep_G(j,4)), rad2deg(kep_G(j,5))];  
end

kep_G(:, 6) = rad2deg(unwrap(kep_G(:, 6)));

%% Filtering
kep_f = movmean(kep_pert, t_span);

%%

figure()

subplot(3, 2, 1)
hold on;
grid on;
plot(t_G, kep_G(:, 1), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert, kep_pert(:, 1), 'DisplayName','Cartesian');
plot(t_pert, kep_f(:, 1), 'b', 'DisplayName', 'Secular');
title('a');
xlabel('time [s]');
ylabel('a [km]');
legend;

subplot(3, 2, 2)
hold on;
grid on;
plot(t_G, kep_G(:, 2), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert, kep_pert(:, 2), 'DisplayName','Cartesian');
plot(t_pert, kep_f(:, 2), 'b', 'DisplayName', 'Secular');
title('e');
xlabel('time [s]');
ylabel('e [-]');
legend;

subplot(3, 2, 3)
hold on;
grid on;
plot(t_G, kep_G(:, 3), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert, kep_pert(:, 3), 'DisplayName','Cartesian');
plot(t_pert, kep_f(:, 3), 'b', 'DisplayName', 'Secular');
title('i');
xlabel('time [s]');
ylabel('i [deg]');
legend;

subplot(3, 2, 4)
hold on;
grid on;
plot(t_G, kep_G(:, 4), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert, kep_pert(:, 4), 'DisplayName','Cartesian');
plot(t_pert, kep_f(:, 4), 'b', 'DisplayName', 'Secular');
title('\Omega');
xlabel('time [s]');
ylabel('\Omega [deg]');
legend;

subplot(3, 2, 5)
hold on;
grid on;
plot(t_G, kep_G(:, 5), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert, kep_pert(:, 5), 'DisplayName','Cartesian');
plot(t_pert, kep_f(:, 5), 'b', 'DisplayName', 'Secular');
title('\omega');
xlabel('time [s]');
ylabel('\omega [deg]');
legend;

subplot(3, 2, 6)
hold on;
grid on;
plot(t_G, kep_G(:, 6), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert, kep_pert(:, 6), 'DisplayName','Cartesian');
plot(t_pert, kep_f(:, 6), 'b', 'DisplayName', 'Secular');
title('\vartheta');
xlabel('time [s]');
ylabel('\vartheta [deg]');
legend;

%% Validation
dOM_theoretical = -(3/2*sqrt(settings.mu)*settings.J2E*settings.RE^2*cos(i)/((1-e^2)^2*a^(7/2)));
dom_theoretical = -(3/2*sqrt(settings.mu)*settings.J2E*settings.RE^2*(5/2*sin(i)^2-2)/((1-e^2)^2*a^(7/2)));

%% Errors

for j = 1:length(t_sample)

    kep_G(j, 3:6) = [deg2rad(kep_G(j,3)), deg2rad(kep_G(j,4)), deg2rad(kep_G(j,5)), deg2rad(kep_G(j,6))];
    kep_pert (j, 3:6) = [deg2rad(kep_pert(j,3)), deg2rad(kep_pert(j,4)), deg2rad(kep_pert(j,5)), deg2rad(kep_pert(j,6))]; 

end

% Relative errors
err = abs(kep_pert-kep_G); 

err(:, 1) = err(:,1)/kep0(1); 
err(:, 3:5) = err(:,3:5)/(2*pi); 
err(:, 6) = err(:,6)./kep_G(:,6); 

figure()

subplot(3,2,1)
semilogy(t_sample/T, err(:, 1))
title('a error');
xlabel('time [T]');
ylabel('a_{car} - a_{gaus} / a_0 [-]');
grid on;

subplot(3,2,2)
semilogy(t_sample/T, err(:,2));
title('e error');
xlabel('time [T]');
ylabel('e_{car} - e_{gauss} [-]');
grid on;

subplot(3,2,3)
semilogy(t_sample/T, err(:,3));
title('i error');
xlabel('time [T]');
ylabel('i_{car} - i_{gauss} / 2\pi [-]');
grid on;

subplot(3,2,4)
semilogy(t_sample/T, err(:,4));
title('\Omega error');
xlabel('time [T]');
ylabel('\Omega_{car} - \Omega_{gauss} / 2\pi [-]');
grid on;

subplot(3,2,5)
semilogy(t_sample/T,err(:,5));
title('\omega error');
xlabel('time [T]');
ylabel('\omega_{car} - \omega_{gauss} / 2\pi [-]');
grid on;

subplot(3,2,6)
semilogy(t_sample/T, err(:,6));
title('\vartheta error');
xlabel('time [T]');
ylabel('\vartheta_{car} - \vartheta_{gauss} / \vartheta_0 [-]');
grid on;

%% Comparison with a real S/C 
% (RBSP A):
% 1 38752U 12046A   22356.68464903  .00041754 -19467-6  20092-2 0  9996
% 2 38752   9.6632 224.8380 6686770  18.4387 357.3584  3.10486121103321

ephemeris = readEphemeris('RBSP_A_1min.txt');

% initial data
a_RBSP = ephemeris(1,1);                 % semi-major axis [km]
e_RBSP = ephemeris(1,2);                 % eccentricity [-]
i_RBSP = deg2rad(ephemeris(1,3));        % inclination [rad]
OM_RBSP = deg2rad(ephemeris(1,4));       % RAAN [rad]
om_RBSP = deg2rad(ephemeris(1,5));       % argument of pericentre [rad]
theta_RBSP = deg2rad(ephemeris(1,6));    % true anomaly [rad]

kep0_RBSP = [a_RBSP e_RBSP i_RBSP OM_RBSP om_RBSP theta_RBSP];

drag_RBSP.c_d = 2.1;                     % drag coefficient [-]
drag_RBSP.Area_mass = 6.64/591.6;        % ratio between reference area and mass [m^2/kg]

T_RBSP = 2*pi*sqrt(a_RBSP^3/settings.mu);
% t_RBSP = 0:1*60:T_RBSP;
% q = length(t_RBSP);

Time = 10*24*3600; % 10 days
t_RBSP = 0:1*60:Time;

options = odeset('RelTol', 1e-10,'AbsTol',1e-11 );
[t_RBSP, kep_RBSP]=ode113(@kep_pert, t_RBSP, kep0_RBSP, options, settings, drag_RBSP);

for j = 1:length(kep_RBSP)

    kep_RBSP(j, 3:6) = [rad2deg(kep_RBSP(j,3)), rad2deg(kep_RBSP(j,4)), ...
            rad2deg(kep_RBSP(j,5)), rad2deg(kep_RBSP(j,6))];  
end

figure()

subplot(3, 2, 1)
hold on;
grid on;
plot(t_RBSP/T_RBSP, kep_RBSP(:, 1), 'DisplayName','modelled a');
plot(t_RBSP/T_RBSP, ephemeris(:, 1), 'DisplayName','real a');
title('a');
xlabel('time [T]');
ylabel('a [km]');
legend;

subplot(3, 2, 2)
hold on;
grid on;
plot(t_RBSP/T_RBSP, kep_RBSP(:, 2), 'DisplayName','modelled e');
plot(t_RBSP/T_RBSP, ephemeris(:, 2), 'DisplayName','real e');
title('e');
xlabel('time [T]');
ylabel('e [-]');
legend;

subplot(3, 2, 3)
hold on;
grid on;
plot(t_RBSP/T_RBSP, kep_RBSP(:, 3), 'DisplayName','modelled i');
plot(t_RBSP/T_RBSP, ephemeris(:, 3), 'DisplayName','real i');
title('i');
xlabel('time [T]');
ylabel('i [deg]');
legend;

subplot(3, 2, 4)
hold on;
grid on;
plot(t_RBSP/T_RBSP, kep_RBSP(:, 4), 'DisplayName','modelled \Omega');
plot(t_RBSP/T_RBSP, ephemeris(:, 4), 'DisplayName','real \Omega');
title('\Omega');
xlabel('time [T]');
ylabel('\Omega [deg]');
legend;

subplot(3, 2, 5)
hold on;
grid on;
plot(t_RBSP/T_RBSP, kep_RBSP(:, 5), 'DisplayName','modelled \omega');
plot(t_RBSP/T_RBSP, ephemeris(:, 5), 'DisplayName','real \omega');
title('\omega');
xlabel('time [T]');
ylabel('\omega [deg]');
legend;

subplot(3, 2, 6)
hold on;
grid on;
plot(t_RBSP/T_RBSP, kep_RBSP(:, 6), 'DisplayName','modelled \vartheta');
plot(t_RBSP/T_RBSP, ephemeris(:, 6), 'DisplayName','real \vartheta');
title('\vartheta');
xlabel('time [T]');
ylabel('\vartheta [deg]');
legend;
