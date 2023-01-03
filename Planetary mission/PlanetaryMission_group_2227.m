%% Orbital Mechanics Project - Planetary mission

%{
    Specifications:
    - semi-major axis = 1.9539*1e4 km
    - eccentricity = 0.6607
    - inclination = 10.4278 deg
    - GT ratio = 13:4
    - Perturbations: J2 and DRAG
    - Aerodynamic drag coefficient = 2.1
    - Area/mass ratio = 0.0042 m^2/kg
    
    Contributors: 
    - Virginia Di Biagio Missaglia
    - Nicolò Galletta
    - Giuseppe Brentino
    - Roberto Pistone Nascone

%}

clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');


%data
keplerian.a = 1.9539*1e4;             % semi-major axis [km]
keplerian.e = 0.6607;                 % eccentricity [-]
keplerian.i = deg2rad(10.4278);       % inclination [rad]
keplerian.OM = deg2rad(60);           % RAAN [rad]
keplerian.om = deg2rad(30);           % argument of pericentre [rad]
keplerian.theta = 0;                  % true anomaly [rad]

drag.c_d = 2.1;             % drag coefficient [-]
drag.Area_mass = 0.0042;    % ratio between reference area and mass [m^2/kg]

settings.mu = astroConstants(13);           % Earth's gravitational constant [km^3/s^2]
settings.J2E = astroConstants(9);           % Earth oblateness parameter [-]
settings.RE = astroConstants(23);           % Earth radius [km]
settings.w_E =2*pi/(23*3600+56*60+4);       % Earth's angular velocity [rad/s]
settings.perturbations = false;

GT_ratio = 13/4;                            % satellite's revs wrt Earth's revs [-]

T = 2*pi*sqrt(keplerian.a^3/settings.mu);   % Period of 1 orbit [s]
T_24 = 23*3600 + 56*60 + 4;                 % 1 day [s] 
T_10 = 10*T_24;                             % 10 days [s]
T_year = 365.25*T_24;                       % 1 year

t0 = 0;
theta_G0 = 0;
t_span = 300;

%% Nominal GT for different orbit periods
[lon_1, lat_1] = GT_data(T, keplerian, settings, drag, theta_G0, t0);
[lon_24, lat_24] = GT_data(T_24, keplerian, settings, drag, theta_G0, t0);
[lon_10, lat_10] = GT_data(T_10, keplerian, settings, drag, theta_G0, t0);

%% Modified GT 
T_mod = (2*pi/settings.w_E)/GT_ratio;
a_mod = ( settings.mu*(T_mod/(2*pi))^2 )^(1/3);

keplerian_mod = keplerian;
keplerian_mod.a = a_mod;

[lon_m, lat_m] = GT_data(T_24*4, keplerian_mod, settings, drag, theta_G0, t0);
[lon_d, lat_d] = GT_data(T_24*4, keplerian, settings, drag, theta_G0, t0);

%% Perturbed GT
settings.perturbations = true;
[lon_pert, lat_pert] = GT_data(T_24, keplerian, settings, drag, theta_G0, t0);

%% Perturbed and Modified GT
[lon_pert_mod, lat_pert_mod] = GT_data(T_mod*13, keplerian_mod, settings, drag, theta_G0, t0);

settings.perturbations = false;
[lon_mod, lat_mod] = GT_data(T_mod*13, keplerian_mod, settings, drag, theta_G0, t0);

%% GT plots

plot_style;

groundTrackPlot
title('1 orbit')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
plot(lon_1, lat_1, 'r');
plot(lon_1(1), lat_1(1), 'o', 'Color', '#EDB120', 'MarkerSize', 12);
plot(lon_1(end), lat_1(end), '*', 'Color', '#EDB120', 'MarkerSize', 12);
legend('groundtrack', 'initial position', 'final position', 'Location','northeast')

groundTrackPlot
title('24 hours orbit')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
plot(lon_24, lat_24, 'r');
plot(lon_24(1), lat_24(1), 'o', 'Color', '#EDB120','MarkerSize', 12);
plot(lon_24(end), lat_24(end), '*', 'Color', '#EDB120', 'MarkerSize', 12);
legend('groundtrack', 'initial position', 'final position', 'Location','northeast')

groundTrackPlot
title('10 days orbit')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
plot(lon_10, lat_10, 'r');
plot(lon_10(1), lat_10(1), 'o', 'Color', '#EDB120', 'MarkerSize', 12);
plot(lon_10(end), lat_10(end), '*', 'Color', '#EDB120', 'MarkerSize', 12);
legend('groundtrack', 'initial position', 'final position', 'Location','northeast')

groundTrackPlot
title('Nominal and Modified ground track')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
plot(lon_m, lat_m, 'g');
plot(lon_m(1), lat_m(1), 'oc', 'MarkerSize', 12);
plot(lon_m(end), lat_m(end), '*c', 'MarkerSize', 12);
plot(lon_d, lat_d, 'r');
plot(lon_d(1), lat_d(1), 'oy', 'MarkerSize', 12);
plot(lon_d(end), lat_d(end), '*y', 'MarkerSize', 12);
legend('modified ground track', 'modified orbit initial position', ...
    'modified orbit final position', 'nominal ground track', ...
    'nominal orbit initial position', 'nominal orbit final position', 'Location','northeast')

groundTrackPlot
title('24 hour orbit')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
plot(lon_pert, lat_pert, 'y');
plot(lon_pert(1), lat_pert(1), 'o', 'Color', '#4DBEEE', 'MarkerSize', 12);
plot(lon_pert(end), lat_pert(end), '*', 'Color', '#4DBEEE', 'MarkerSize', 12);
plot(lon_24, lat_24, 'r');
plot(lon_24(1), lat_24(1), 'or', 'MarkerSize', 12);
plot(lon_24(end), lat_24(end), '*r', 'MarkerSize', 12);
legend('perturbed gt', 'initial perturbed position', ...
    'final perturbed position', 'nominal gt', 'initial nominal position', ...
    'final nominal position', 'Location','northeast')

groundTrackPlot
title('13 orbits')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
plot(lon_pert_mod, lat_pert_mod, 'm');
plot(lon_pert_mod(1), lat_pert_mod(1), 'o', 'Color', '#EDB120', 'MarkerSize', 12);
plot(lon_pert_mod(end), lat_pert_mod(end), '*', 'Color', '#EDB120', 'MarkerSize', 12);
plot(lon_mod, lat_mod, 'g');
plot(lon_mod(1), lat_mod(1), 'or', 'MarkerSize', 12);
plot(lon_mod(end), lat_mod(end), '*r', 'MarkerSize', 12);
legend('perturbed and modified gt', 'initial pert. and mod. position', ...
    'final pert. and mod. position', 'modified gt', 'initial mod. position', ...
    'final mod. position', 'Location','northeast')


%% Orbit propagation with Gauss's Eqns

settings.perturbations = true;
t_sample = 0:T/t_span:T_24*180;

options = odeset('RelTol', 1e-13,'AbsTol',1e-14 );

kep0 = [keplerian.a keplerian.e keplerian.i keplerian.OM keplerian.om keplerian.theta];

tic
[t_G, kep_G] = ode113(@GaussEq_rsw, t_sample, kep0, options, settings, drag);
t_Gauss = toc;

[r0, v0] = kep2car(keplerian.a, keplerian.e, keplerian.i, keplerian.OM, ... 
    keplerian.om, keplerian.theta, settings.mu);
s0 = [r0; v0];

tic
[t_pert, Y_pert] = ode113(@pert_tbp, t_sample, s0, options, settings, drag);
t_Cart = toc;

% figure()
% title('Orbit plot - Cartesian Propagation')
% hold on;
% grid on;
% axis equal;
% xlabel('x [km]')
% ylabel('y [km]')
% zlabel('z [km]')
% plot_terra;
% plot3(Y_pert(:, 1), Y_pert(:, 2), Y_pert(:, 3));

%% Errors

kep_pert = zeros(size(Y_pert,1),6);

for j = 1:size(Y_pert,1)

    r = Y_pert(j, 1:3);
    v = Y_pert(j, 4:6);

    [a_p, e_p, i_p, OM_p, om_p, th_p] = car2kep(r, v, settings.mu);
    kep_pert(j, :) = [a_p, e_p, i_p, OM_p, om_p, th_p]; 

end

kep_pert(:, 4) = wrapTo2Pi(kep_pert(:, 4));
kep_pert(:, 5) = wrapTo2Pi(kep_pert(:, 5));
kep_pert(:, 6) = wrapTo360(kep_pert(:, 6));

kep_G(:, 4) = wrapTo2Pi(kep_G(:, 4));
kep_G(:, 5) = wrapTo2Pi(kep_G(:, 5));
kep_G(:, 6) = wrapTo2Pi(kep_G(:, 6));

% comparson of the two methods

err = abs(kep_pert-kep_G); 

err(:, 1) = err(:,1)/kep0(1); 
err(:, 3:5) = err(:,3:5)/(2*pi); 
err(:, 6) = err(:,6)./kep_G(:,6); 

figure('Units','normalized', 'OuterPosition',[0 0 1 1])

subplot(3,2,1)
semilogy(t_sample/T, err(:, 1))
title('a error');
xlabel('time [T]');
ylabel('$|a_{car}$ - $a_{gauss}|$ / $a_0$ [-]');
grid on;

subplot(3,2,2)
semilogy(t_sample/T, err(:,2));
title('e error');
xlabel('time [T]');
ylabel('$|e_{car}$ - $e_{gauss}|$ [-]');
grid on;

subplot(3,2,3)
semilogy(t_sample/T, err(:,3));
title('i error');
xlabel('time [T]');
ylabel('$|i_{car}$ - $i_{gauss}|$ / $2\pi$ [-]');
grid on;

subplot(3,2,4)
semilogy(t_sample/T, err(:,4));
title('$\Omega$ error');
xlabel('time [T]');
ylabel('$|\Omega_{car}$ - $$\Omega_{gauss}|$$ / $2\pi$ [-]');
grid on;

subplot(3,2,5)
semilogy(t_sample/T,err(:,5));
title('$\omega$ error');
xlabel('time [T]');
ylabel('$$|\omega_{car}$$ - $$\omega_{gauss}|$$ / $2\pi$ [-]');
grid on;

subplot(3,2,6)
semilogy(t_sample/T, err(:,6));
title('$\vartheta$ error');
xlabel('time [T]');
ylabel('$$|\vartheta_{car}$$ - $$\vartheta_{gauss}|$$ / $$\vartheta_{gauss}$$ [-]');
grid on;


%% History of Keplerian elements

for j = 1:length(kep_G)

    kep_G(j, 3:6) = [rad2deg(kep_G(j,3)), rad2deg(kep_G(j,4)), rad2deg(kep_G(j,5)), rad2deg(kep_G(j,6))]; 
    kep_pert(j, 3:6) = [rad2deg(kep_pert(j,3)), rad2deg(kep_pert(j,4)), rad2deg(kep_pert(j,5)), rad2deg(kep_pert(j,6))];

end


%% Filtering
kep_f = movmean(kep_pert, t_span);

%% Figures

figure('Units','normalized', 'OuterPosition',[0 0 1 1])

subplot(3, 2, 1)
hold on;
grid on;
plot(t_G/T, kep_G(:, 1), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert/T, kep_pert(:, 1), 'DisplayName','Cartesian');
plot(t_pert/T, kep_f(:, 1), 'b', 'DisplayName', 'Secular');
title('a');
xlabel('time [T]');
ylabel('a [km]');

subplot(3, 2, 2)
hold on;
grid on;
plot(t_G/T, kep_G(:, 2), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert/T, kep_pert(:, 2), 'DisplayName','Cartesian');
plot(t_pert/T, kep_f(:, 2), 'b', 'DisplayName', 'Secular');
title('e');
xlabel('time [T]');
ylabel('e [-]');

subplot(3, 2, 3)
hold on;
grid on;
plot(t_G/T, kep_G(:, 3), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert/T, kep_pert(:, 3), 'DisplayName','Cartesian');
plot(t_pert/T, kep_f(:, 3), 'b', 'DisplayName', 'Secular');
title('i');
xlabel('time [T]');
ylabel('i [deg]');

subplot(3, 2, 4)
hold on;
grid on;
plot(t_G(1:300*5)/T, kep_G(1:300*5, 4), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert(1:300*5)/T, kep_pert(1:300*5, 4), 'DisplayName','Cartesian');
plot(t_pert(1:300*5)/T, kep_f(1:300*5, 4), 'b', 'DisplayName', 'Secular');
title('$\Omega$');
xlabel('time [T]');
ylabel('$\Omega$ [deg]');
legend;

subplot(3, 2, 5)
hold on;
grid on;
plot(t_G(1:300*5)/T, kep_G(1:300*5, 5), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert(1:300*5)/T, kep_pert(1:300*5, 5), 'DisplayName','Cartesian');
plot(t_pert(1:300*5)/T, kep_f(1:300*5, 5), 'b', 'DisplayName', 'Secular');
title('$\omega$');
xlabel('time [T]');
ylabel('$\omega$ [deg]');

subplot(3, 2, 6)
hold on;
grid on;
plot(t_G(1:300*5)/T, kep_G(1:300*5, 6), 'DisplayName','Gauss Equation(RSW)');
plot(t_pert(1:300*5)/T, kep_pert(1:300*5, 6), 'DisplayName','Cartesian');
plot(t_pert(1:300*5)/T, kep_f(1:300*5, 6), 'b', 'DisplayName', 'Secular');
title('$\vartheta$');
xlabel('time [T]');
ylabel('$\vartheta$ [deg]');

%% Validation

dOM_theoretical = -(3/2*sqrt(settings.mu)*settings.J2E*settings.RE^2*...
    cos(keplerian.i)/((1-keplerian.e^2)^2*keplerian.a^(7/2)));
dom_theoretical = -(3/2*sqrt(settings.mu)*settings.J2E*settings.RE^2*...
    (5/2*sin(keplerian.i)^2-2)/((1-keplerian.e^2)^2*keplerian.a^(7/2)));

%% Comparison with a real S/C 
% (RBSP A):
% 1 38752U 12046A   22356.68464903  .00041754 -19467-6  20092-2 0  9996
% 2 38752   9.6632 224.8380 6686770  18.4387 357.3584  3.10486121103321
ephemeris = readEphemeris('RBSP_A_1min_10d');

% initial data
a_SC = ephemeris(1,1);                 % semi-major axis [km]
e_SC = ephemeris(1,2);                 % eccentricity [-]
i_SC = deg2rad(ephemeris(1,3));        % inclination [rad]
OM_SC = deg2rad(ephemeris(1,4));       % RAAN [rad]
om_SC = deg2rad(ephemeris(1,5));       % argument of pericentre [rad]
theta_SC = deg2rad(ephemeris(1,6));    % true anomaly [rad]

kep0_SC = [a_SC e_SC i_SC OM_SC om_SC theta_SC];

drag_SC.c_d = 2.1;                     % drag coefficient [-]
drag_SC.Area_mass = (1.8*1.8+4*0.9*0.9) /647.6; % ratio between reference area and mass [m^2/kg]

T_SC = 2*pi*sqrt(a_SC^3/settings.mu);

Time = 100*T_SC;
t_span_SC = 0:T_SC/300:Time;
t_SC_ephem = 0:60:10*24*3600;

options = odeset('RelTol', 1e-13,'AbsTol',1e-14 );
settings.perturbations = true;
[t_SC_kep, kep_SC]=ode113(@GaussEq_rsw, t_span_SC, kep0_SC, options, settings, drag_SC);

for j=1:length(t_SC_kep)
    kep_SC(j, 3:6) = [rad2deg(kep_SC(j,3)), rad2deg(kep_SC(j,4)), ...
        rad2deg(kep_SC(j,5)), rad2deg(wrapTo2Pi(kep_SC(j,6)))];
end

% plots
figure()

subplot(3, 2, 1)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 1), 'DisplayName','modelled a');
plot(t_SC_ephem/T_SC, ephemeris(:, 1), 'DisplayName','real a');
title('a');
xlabel('time [T]');
ylabel('a [km]');
legend;

subplot(3, 2, 2)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 2), 'DisplayName','modelled e');
plot(t_SC_ephem/T_SC, ephemeris(:, 2), 'DisplayName','real e');
title('e');
xlabel('time [T]');
ylabel('e [-]');
legend;

subplot(3, 2, 3)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 3), 'DisplayName','modelled i');
plot(t_SC_ephem/T_SC, ephemeris(:, 3), 'DisplayName','real i');
title('i');
xlabel('time [T]');
ylabel('i [deg]');
legend;

subplot(3, 2, 4)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 4), 'DisplayName','modelled \Omega');
plot(t_SC_ephem/T_SC, ephemeris(:, 4), 'DisplayName','real \Omega');
title('\Omega');
xlabel('time [T]');
ylabel('\Omega [deg]');
legend;

subplot(3, 2, 5)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 5), 'DisplayName','modelled \omega');
plot(t_SC_ephem/T_SC, ephemeris(:, 5), 'DisplayName','real \omega');
title('\omega');
xlabel('time [T]');
ylabel('\omega [deg]');
legend;

subplot(3, 2, 6)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 6), 'DisplayName','modelled \vartheta');
plot(t_SC_ephem/T_SC, ephemeris(:, 6), 'DisplayName','real \vartheta');
title('\vartheta');
xlabel('time [T]');
ylabel('\vartheta [deg]');
legend;