%% Comparison with a real S/C 
% (RBSP A):
% 1 38752U 12046A   22356.68464903  .00041754 -19467-6  20092-2 0  9996
% 2 38752   9.6632 224.8380 6686770  18.4387 357.3584  3.10486121103321

clearvars; close all; clc

addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');

settings.mu = astroConstants(13);           % Earth's gravitational constant [km^3/s^2]
settings.J2E = astroConstants(9);           % Earth oblateness parameter [-]
settings.RE = astroConstants(23);           % Earth radius [km]
settings.w_E = deg2rad(15.04)/3600;         % Earth's angular velocity [rad/s]

ephemeris = readEphemeris('RBSP_A_10day_1min.txt');

% initial data)
a_SC = ephemeris(1,1);                 % semi-major axis [km]
e_SC = ephemeris(1,2);                 % eccentricity [-]
i_SC = deg2rad(ephemeris(1,3));        % inclination [rad]
OM_SC = deg2rad(ephemeris(1,4));       % RAAN [rad]
om_SC = deg2rad(ephemeris(1,5));       % argument of pericentre [rad]
theta_SC = deg2rad(ephemeris(1,6));    % true anomaly [rad]

kep0_SC = [a_SC e_SC i_SC OM_SC om_SC theta_SC];
[un, dos] = kep2car(a_SC, e_SC, i_SC, OM_SC, om_SC, theta_SC);
s0_SC = [un; dos];

drag_SC.c_d = 2.1;                     % drag coefficient [-]

% ratio between reference area and mass [m^2/kg]
% drag_SC.Area_mass = 6.64 /591.6; % RBSP_A 
drag_SC.Area_mass = 0.3142/28.7; % vesselSat
% drag_SC.Area_mass = 0.6636/104; % ov1
% drag_SC.Area_mass = 0.817/78.5; % ops1427

T_SC = 2*pi*sqrt(a_SC^3/settings.mu);

Time = 10*24*3600; % 10 days
t_SC_ephem = 0:1*60:Time;

options = odeset('RelTol', 1e-10,'AbsTol',1e-11 );
settings.perturbations = true;
[t_SC_car, Y_SC]=ode113(@pert_tbp, [0 Time], s0_SC, options, settings, drag_SC);

kep_SC_car = zeros(size(Y_SC));
for j= 1:length(t_SC_car)
    [un, dos, tres, quatros, cinq, sei] = car2kep(Y_SC(j,1:3), Y_SC(j, 4:6), settings.mu);
    kep_SC_car(j, :) = [un, dos, tres, quatros, cinq, sei];
    kep_SC_car(j, 3:6) = [rad2deg(kep_SC_car(j, 3)), rad2deg(kep_SC_car(j, 4)), ...
        rad2deg(kep_SC_car(j, 5)), rad2deg(kep_SC_car(j, 6))];
end

[t_SC_kep, kep_SC]=ode113(@kep_pert, [0 Time], kep0_SC, options, settings, drag_SC);

%% Accelerations
settings.ref_sys = 'RSW';

a_p=zeros(3,length(t_SC_kep));
a_J2=zeros(3,length(t_SC_kep));
a_drag=zeros(3,length(t_SC_kep));
r = zeros(length(t_SC_kep), 1);

for j=1:length(t_SC_kep)
    r(j) = norm(kep2car(kep_SC(j,1), kep_SC(j,2), kep_SC(j,3), kep_SC(j,4), kep_SC(j,5), kep_SC(j,6)));
    [a_p(:,j), a_J2(:,j), a_drag(:,j)] = a_pert(settings, kep_SC(j,:), drag_SC);
end

r = r -settings.RE;

figure()
plot(t_SC_kep, r)

figure()
plot(t_SC_kep, a_J2(1, :), t_SC_kep, a_J2(2, :), t_SC_kep, a_J2(3, :))
legend('a_{J2} r', 'a_{J2} s', 'a_{J2} w')

figure()
plot(r, a_p(1, :), r, a_p(2, :), r, a_p(3, :))
legend('a_p r', 'a_p s', 'a_p w')

figure()
plot(r, a_J2(1, :), r, a_J2(2, :), r, a_J2(3, :))
legend('a_{J2} r', 'a_{J2} s', 'a_{J2} w')

figure()
plot(r, a_drag(1, :), r, a_drag(2, :), r, a_drag(3, :))
legend('a_{drag} r', 'a_{drag} s', 'a_{drag} w')


%% Subplot

for j = 1:length(t_SC_kep)

    kep_SC(j, 3:6) = [rad2deg(wrapToPi(kep_SC(j,3))), rad2deg(wrapTo2Pi(kep_SC(j,4))), ...
            rad2deg(wrapTo2Pi(kep_SC(j,5))), rad2deg(wrapTo2Pi(kep_SC(j,6)))];  
end

figure()

subplot(3, 2, 1)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 1), 'DisplayName','modelled a');
plot(t_SC_car/T_SC, kep_SC_car(:, 1), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 1), 'DisplayName','real a');
title('a');
xlabel('time [T]');
ylabel('a [km]');
legend;

subplot(3, 2, 2)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 2), 'DisplayName','modelled e');
plot(t_SC_car/T_SC, kep_SC_car(:, 2), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 2), 'DisplayName','real e');
title('e');
xlabel('time [T]');
ylabel('e [-]');
legend;

subplot(3, 2, 3)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 3), 'DisplayName','modelled i');
plot(t_SC_car/T_SC, kep_SC_car(:, 3), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 3), 'DisplayName','real i');
title('i');
xlabel('time [T]');
ylabel('i [deg]');
legend;

subplot(3, 2, 4)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 4), 'DisplayName','modelled \Omega');
plot(t_SC_car/T_SC, kep_SC_car(:, 4), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 4), 'DisplayName','real \Omega');
title('\Omega');
xlabel('time [T]');
ylabel('\Omega [deg]');
legend;

subplot(3, 2, 5)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 5), 'DisplayName','modelled \omega');
plot(t_SC_car/T_SC, kep_SC_car(:, 5), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 5), 'DisplayName','real \omega');
title('\omega');
xlabel('time [T]');
ylabel('\omega [deg]');
legend;

subplot(3, 2, 6)
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 6), 'DisplayName','modelled \vartheta');
plot(t_SC_car/T_SC, kep_SC_car(:, 6), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 6), 'DisplayName','real \vartheta');
title('\vartheta');
xlabel('time [T]');
ylabel('\vartheta [deg]');
legend;

%% plot singoli
figure()
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 1), 'DisplayName','modelled a');
plot(t_SC_car/T_SC, kep_SC_car(:, 1), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 1), 'DisplayName','real a');
title('a');
xlabel('time [T]');
ylabel('a [km]');
legend;

figure()
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 2), 'DisplayName','modelled e');
plot(t_SC_car/T_SC, kep_SC_car(:, 2), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 2), 'DisplayName','real e');
title('e');
xlabel('time [T]');
ylabel('e [-]');
legend;

figure()
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 3), 'DisplayName','modelled i');
plot(t_SC_car/T_SC, kep_SC_car(:, 3), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 3), 'DisplayName','real i');
title('i');
xlabel('time [T]');
ylabel('i [deg]');
legend;

figure()
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 4), 'DisplayName','modelled \Omega');
plot(t_SC_car/T_SC, kep_SC_car(:, 4), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 4), 'DisplayName','real \Omega');
title('\Omega');
xlabel('time [T]');
ylabel('\Omega [deg]');
legend;

figure()
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 5), 'DisplayName','modelled \omega');
plot(t_SC_car/T_SC, kep_SC_car(:, 5), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 5), 'DisplayName','real \omega');
title('\omega');
xlabel('time [T]');
ylabel('\omega [deg]');
legend;

figure()
hold on;
grid on;
plot(t_SC_kep/T_SC, kep_SC(:, 6), 'DisplayName','modelled \vartheta');
plot(t_SC_car/T_SC, kep_SC_car(:, 6), 'DisplayName','modelled a car');
plot(t_SC_ephem/T_SC, ephemeris(:, 6), 'DisplayName','real \vartheta');
title('\vartheta');
xlabel('time [T]');
ylabel('\vartheta [deg]');
legend;