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
kep0 = [a e deg2rad(i) deg2rad(OM) deg2rad(om) deg2rad(theta)];

GT_ratio = 13/4;            % satellite's revs wrt Earth's revs [-]
drag.c_d = 2.1;                  % drag coefficient [-]
drag.Area_mass = 0.0042;         % ratio between reference area and mass [m^2/kg]

settings.mu = astroConstants(13);               % Earth's gravitational constant [km^3/s^2]
settings.J2E = astroConstants(9);               % Earth oblateness parameter [-]
settings.RE = astroConstants(23);               % Earth radius [km]
settings.w_E = deg2rad(15.04)/3600;              % Earth's angular velocity [rad/s]


T = 2*pi*sqrt(a^3/settings.mu);             % Period of 1 orbit [s]
T_24 = 24*3600;                             % 1 day [s] 
T_10 = 10*T_24;                             % 10 days [s]


%% nominal GT 

%T/T_24;
[r0, v0] = kep2car(a, e, i, OM, om, theta, settings.mu);
s0 = [r0; v0];                              % initial state vector

settings.perturbations = false;
options = odeset('RelTol', 1e-10,'AbsTol',1e-11 );
[t, Y] = ode113(@pert_tbp,[0 T], s0, options, settings, drag);


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

% TOGLIERE I MARGINI BIANCHI DAL PLOT DELLA GT
groundTrackPlot
% xlim([-180 180]);
% ylim([-90 90]);
plot(lon, lat, 'r','LineWidth',2);
plot(lon(1), lat(1), 'o', 'LineWidth', 5);
plot(lon(end), lat(end), '*', 'LineWidth', 5);


%% Modified GT

T_mod = (2*pi/settings.w_E)/GT_ratio;
a_mod = ( settings.mu*(T_mod/(2*pi))^2 )^(1/3);

[r_mod, v_mod] = kep2car(a_mod, e, i, OM, om, theta, settings.mu);
s_mod = [r_mod; v_mod];                              % initial modified state vector

settings.perturbations = false;
[t1, Y1] = ode113(@pert_tbp,[0 T_mod], s_mod, options, settings, drag);

[alpha_mod, delta_mod, lon_mod, lat_mod] = groundTrack (Y1, theta_G0, t1, settings.w_E, t0);

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
tic

settings.perturbations = true;
[t_pert, Y_pert] = ode113(@pert_tbp,[0 T], s0, options, settings, drag);

figure()
hold on;
axis equal;
plot_terra;
plot3(Y(:, 1), Y(:, 2), Y(:, 3));
plot3(Y_pert(:, 1), Y_pert(:, 2), Y_pert(:, 3));

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


% TOGLIERE I MARGINI BIANCHI DAL PLOT DELLA GT
groundTrackPlot
% xlim([-180 180]);
% ylim([-90 90]);
plot(lon_p, lat_p, 'y','LineWidth',2);
plot(lon_p(1), lat_p(1), 'o', 'LineWidth', 5);
plot(lon_p(end), lat_p(end), '*', 'LineWidth', 5);

% IN ITALIANO: considerate le nobili origini della mamma di MATLAB,
% osserviamo che la GT plottata ESATTAMENTE per il tempo T di un'orbita
% restituisce una BELLISSIMA LINEA DRITTA tra penultimo e ultimo punto (non
% c'è nulla in mezzo)

toc
%% Perturbated and Modified GT
tic

settings.perturbations = true;   
[t_m_pert, Y_m_pert] = ode113(@pert_tbp,[0 T_mod], s_mod, options, settings, drag);

[alpha_p_mod, delta_p_mod, lon_p_mod, lat_p_mod] = groundTrack (Y_m_pert, ...
    theta_G0, t_m_pert, settings.w_E, t0);

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

% TOGLIERE I MARGINI BIANCHI DAL PLOT DELLA GT
groundTrackPlot
% xlim([-180 180]);
% ylim([-90 90]);
plot(lon_p_mod, lat_p_mod, 'y','LineWidth',2);
plot(lon_p_mod(1), lat_p_mod(1), 'o', 'LineWidth', 5);
plot(lon_p_mod(end), lat_p_mod(end), '*', 'LineWidth', 5);


toc


%% Orbit propagation with Gauss's Eqns

[t_G, Y_G] = ode113(@kep_pert, [0 T], kep0, options, settings, drag);

for j = 1:size(Y_G, 1)
    [Y_G(j, 1:3), Y_G(j, 4:6)] = kep2car(Y_G(j, 1), Y_G(j, 2), Y_G(j, 3), Y_G(j, 4), Y_G(j, 5), Y_G(j, 6));

end

figure()
hold on;
axis equal;
plot_terra;
plot3(Y_pert(:, 1), Y_pert(:, 2), Y_pert(:, 3));
plot3(Y_G(:, 1), Y_G(:, 2), Y_G(:, 3));
legend('Orbit car', 'Orbit kep');

% %% History of Keplerian elements
% 
for j = 1:size(Y_pert,1)

    r = Y_pert(j, 1:3);
    v = Y_pert(j, 4:6);

    [a, e, i, OM, om, th] = car2kep(r, v, settings.mu);
    kep(j, :) = [a, e, rad2deg(i), rad2deg(OM), rad2deg(om), rad2deg(th)];

end

%kep(:, 6) = rad2deg(unwrap(kep(:, 6)));
% 
% for j = 1:length(Y_G)
% 
%     Y_G(j, 3:6) = [rad2deg(Y_G(j,3)), rad2deg(Y_G(j,4)), rad2deg(Y_G(j,5)), rad2deg(Y_G(j,6))];  
% 
% end
% 
% 
% t_car = t_pert;
% t_RSW = t_G;
% 
% figure()
% 
% subplot(3, 2, 1)
% hold on;
% grid on;
% plot(t_RSW, Y_G(:, 1), 'DisplayName','Gauss Equation(RSW)');
% plot(t_car, kep(:, 1), 'DisplayName','Cartesian');
% %plot(t, kep_f(:, 1), 'DisplayName', 'Secular');
% title('a');
% xlabel('time [T]');
% ylabel('a [km]');
% legend;
% 
% subplot(3, 2, 2)
% hold on;
% grid on;
% plot(t_RSW, Y_G(:, 2), 'DisplayName','Gauss Equation(TNH)');
% plot(t_car, kep(:, 2), 'DisplayName','Cartesian');
% %plot(t, kep_f(:, 2), 'DisplayName', 'Secular');
% title('e');
% xlabel('time [T]');
% ylabel('e [-]');
% legend;
% 
% subplot(3, 2, 3)
% hold on;
% grid on;
% plot(t_RSW, Y_G(:, 3), 'DisplayName','Gauss Equation(TNH)');
% plot(t_car, kep(:, 3), 'DisplayName','Cartesian');
% %plot(t, kep_f(:, 3), 'DisplayName', 'Secular');
% title('i');
% xlabel('time [T]');
% ylabel('i [deg]');
% legend;
% 
% subplot(3, 2, 4)
% hold on;
% grid on;
% plot(t_RSW, Y_G(:, 4), 'DisplayName','Gauss Equation(TNH)');
% plot(t_car, kep(:, 4), 'DisplayName','Cartesian');
% %plot(t, kep_f(:, 4), 'DisplayName', 'Secular');
% title('\Omega');
% xlabel('time [T]');
% ylabel('\Omega [deg]');
% legend;
% 
% subplot(3, 2, 5)
% hold on;
% grid on;
% plot(t_RSW, Y_G(:, 5), 'DisplayName','Gauss Equation(TNH)');
% plot(t_car, kep(:, 5), 'DisplayName','Cartesian');
% %plot(t, kep_f(:, 5), 'DisplayName', 'Secular');
% title('\omega');
% xlabel('time [T]');
% ylabel('\omega [deg]');
% legend;
% 
% subplot(3, 2, 6)
% hold on;
% grid on;
% plot(t_RSW, Y_G(:, 6), 'DisplayName','Gauss Equation(TNH)');
% plot(t_car, kep(:, 6), 'DisplayName','Cartesian');
% %plot(t, kep_f(:, 6), 'DisplayName', 'Secular');
% title('\vartheta');
% xlabel('time [T]');
% ylabel('\vartheta [deg]');
% legend;




