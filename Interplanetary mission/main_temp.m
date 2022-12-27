%% Orbital Mechanics Project - Interplanetary mission

clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('../shared functions/planet3D/');
addpath('./functions/');

% define gravitational constants
mu_Sun = astroConstants(4);
mu_E = astroConstants(13);
mu_Sat = astroConstants(16);

% define indexes for uplanet and ephNEO functions

i_flyby = 6;    % Saturn
i_NEO = 21;     % NEO
i_dep = 3; % Earth

%% Saturn-NEO Leg

% define MJD2000 dates for the Saturn-NEO heliocentric leg, with almost
% every date available and a coarse grid 

S_N.departure.earliest = [2031, 8, 2, 0, 0, 0]; 
S_N.departure.latest   = [2064, 1, 28, 23, 59, 59];

S_N.arrival.earliest = [2032, 8, 2, 0, 0, 0];
S_N.arrival.latest   = [2065, 1, 28, 23, 59, 59];

delta_t = days(90); % days separating to dates in the search grid

[dates.flyby] = timeWindowMj2000(S_N.departure, delta_t);
[dates.arrival] = timeWindowMj2000(S_N.arrival, delta_t);

% find the required deltaV for each possible transfer
m = length(dates.flyby);
n = length(dates.arrival);

S_N.deltaV = zeros(m, n);

S_N.deltaV_dep = cell(m, n);
S_N.VI = cell(m, n);

S_N.Tpar= zeros(m, n);
for i = 1:m
    for j = 1:n
        [S_N.deltaV_arr(i, j), S_N.deltaV_dep{i, j}, S_N.VI{i,j}, S_N.Tpar(i, j)] = ...
            SNEO_transfer (dates.flyby(i), dates.arrival(j), mu_Sun, i_flyby, i_NEO);
         S_N.deltaV(i,j) = S_N.deltaV_arr(i,j) + S_N.deltaV_dep{i,j};
    end
end


%% Cheapest dates departure/arrival to NEO  

dates.flyby_logic = false(length(dates.flyby), length(dates.arrival));

i = 1;
 for k = 1:size(S_N.deltaV_arr, 1)
    for j = 1:size(S_N.deltaV_arr, 2)
        temp.time_arr = datetime(mjd20002date(dates.arrival(j)));
         temp.time_dep = datetime(mjd20002date(dates.flyby(k)));

        if S_N.deltaV(k, j) < 40 && S_N.Tpar(k, j) < ( seconds(diff([temp.time_dep; temp.time_arr] ) ) ) % S_N.deltaV prima era S_N.deltaV_arr
            
            dates.flyby_logic(k, j) = true;

            S_N.data{i} = [dates.flyby(k); dates.arrival(j);...
                S_N.deltaV(k,j); S_N.VI{k,j}'; S_N.deltaV_dep{k,j} ];

            i = i+1;

        end
    end
 end


%% PorkChop plot Saturn-NEO

figure()
hold on;
contour(dates.flyby, dates.arrival, S_N.deltaV_arr', 20:2:40, 'LineWidth',2); % S_N.deltaV prima era S_N.deltaV_arr
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([20, 40]);

%% Save best flyby dates and correspondent deltaV_arr
temp.best_fb_time = find(~all(dates.flyby_logic==0,2)); % find the rows with min delta v (fb date)

dates.flyby = dates.flyby(temp.best_fb_time);

S_N.deltaV_arr = S_N.deltaV_arr(temp.best_fb_time,:); % extract delta_v rows of best fb date probably useless

%% Earth-Saturn Leg

E_S.departure.earliest = [2030, 8, 2, 0, 0, 0]; 

dates.departure = linspace(min(dates.flyby(1)-730, date2mjd2000(E_S.departure.earliest) ),...
     dates.flyby(end)-180, 70);

m = length(dates.departure);
n = length(dates.flyby);

% indeces i, m for dates.departure elements
% indeces j, n for dates.flyby elements

E_S.deltaV_dep = zeros(m, n);
E_S.VF = cell(m, n);
E_S.VI = cell(m, n);
E_S.T_par = zeros(m, n);

for i = 1:m
    for j = 1:n
        [E_S.deltaV_dep(i, j), E_S.VF{i, j}, E_S.deltaVI(i, j), E_S.T_par(i, j)] =...
            ES_transfer (dates.departure(i), dates.flyby(j), mu_Sun, i_dep, i_flyby);
      
    end
end


%% PorkChop plot Saturn-Earth

figure()
hold on;
contour(dates.departure, dates.flyby, E_S.deltaV_dep', 20:2:30, 'LineWidth',2);
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([20, 30]);

%% Cheapest departure/arrival to Saturn dates 

i = 1;
for k = 1:size(E_S.deltaV_dep, 1)
    for j = 1:size(E_S.deltaV_dep, 2)
 
        temp.time_fb = datetime(mjd20002date(dates.flyby(j)));
        temp.time_dep = datetime(mjd20002date(dates.departure(k)));

        if E_S.deltaV_dep(k, j) < 25 && E_S.T_par(k, j) < ( seconds(diff([temp.time_dep; temp.time_fb] ) ) )
            E_S.data{i} = [dates.departure(k); dates.flyby(j); E_S.deltaV_dep(k,j);...
                E_S.deltaVI(k,j); E_S.VF{k,j}'];
            i = i+1;

        end
    end
end

%% Find best date
k = 1;
for i = 1 : length(S_N.data)
    for j = 1 : length(E_S.data)
        if S_N.data{i}(1) == E_S.data{j}(2)
            temp.deltaV = E_S.data{j}(4) + norm(E_S.data{j}(5:7)-S_N.data{i}(4:6))...
                + S_N.data{i}(7);
            mission_data(:,k) =[E_S.data{j}(1); E_S.data{j}(2); S_N.data{i}(2); temp.deltaV];
            k = k+1;
        end
    end
end

[~,index] = min(mission_data(4,:));

%% Optimization

guess = [mission_data(1,index);mission_data(2,index);mission_data(3,index)];

lb = guess - 180; 
ub = guess + 180; 
 
options = optimoptions("fmincon", 'Algorithm','sqp');
[u, DeltaV_opt, ~, ~] = fmincon(@(u)obj_Fcn(u,i_dep,i_flyby,i_NEO),...
    guess, [], [], [], [], lb, ub,@(u) nonlincon(u,i_dep,i_flyby,i_NEO), options);

dates.departure = u(1);
dates.flyby  = u(2);
dates.arrival = u(3);

%% Heliocentric trajectory plot

[plot.kep_sat, ~] = uplanet(dates.flyby, i_flyby);
[plot.r_sat, ~] = kep2car(plot.kep_sat(1),plot.kep_sat(2),plot.kep_sat(3),...
    plot.kep_sat(4),plot.kep_sat(5),plot.kep_sat(6),mu_Sun);

[plot.kep_earth, ~] = uplanet(dates.departure, i_dep);
[plot.r_earth, ~] = kep2car(plot.kep_earth(1),plot.kep_earth(2),plot.kep_earth(3),...
    plot.kep_earth(4),plot.kep_earth(5),plot.kep_earth(6),mu_Sun);

[plot.kep_NEO, ~] = ephNEO(dates.arrival, i_NEO);
[plot.r_NEO, ~] = kep2car(plot.kep_NEO(1),plot.kep_NEO(2),plot.kep_NEO(3),...
    plot.kep_NEO(4),plot.kep_NEO(5),plot.kep_NEO(6),mu_Sun);

dates.dep = datetime(mjd20002date(dates.departure));
dates.fb = datetime(mjd20002date(dates.flyby));
dates.arr = datetime(mjd20002date(dates.arrival));

E_S.TOF = seconds(diff([dates.dep; dates.fb]));
S_N.TOF = seconds(diff([dates.fb; dates.arr]));


[~,~,~,~,plot.VI_e,~,~,~] = lambertMR(plot.r_earth, plot.r_sat,E_S.TOF, mu_Sun, 0, 0, 2);
[~,~,~,~,plot.VI_s,~,~,~] = lambertMR(plot.r_sat, plot.r_NEO,S_N.TOF, mu_Sun, 0, 0, 2);

plot.s0_e = [plot.r_earth; plot.VI_e'];
plot.s0_s = [plot.r_sat; plot.VI_s'];

settings.perturbations = false;
settings.mu = mu_Sun;
options = odeset('RelTol', 1e-10,'AbsTol',1e-11);

[t_e, Y_e] = ode113(@pert_tbp, [0 E_S.TOF], plot.s0_e, options, settings);
[t_s, Y_s] = ode113(@pert_tbp, [0 S_N.TOF], plot.s0_s, options, settings);


figure('Name', 'Lambert')
hold on
plot3( Y_e(:, 1), Y_e(:, 2), Y_e(:, 3), '-')
plot3( Y_s(:, 1), Y_s(:, 2), Y_s(:, 3), '-')
plot3(plot.r_earth(1), plot.r_earth(2), plot.r_earth(3), '.', 'MarkerSize', 20);
plot3(plot.r_sat(1), plot.r_sat(2), plot.r_sat(3), '.', 'MarkerSize', 20);
plot3(plot.r_NEO(1), plot.r_NEO(2), plot.r_NEO(3),'.', 'MarkerSize', 20);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Transfer Orbit with Lambert Algorithm')
axis equal;
grid on;

    
%% Saturn Fly-by

% saturn parameters
[kep_saturn,~] = uplanet( dates.flyby, i_flyby);

[r_saturn, v_saturn] = kep2car(kep_saturn(1), kep_saturn(2), kep_saturn(3),...
    kep_saturn(4), kep_saturn(5), kep_saturn(6), mu_Sun);

time_fb = datetime( mjd20002date( dates.flyby ) );

% earth parameters
[kep_dep, ~] = uplanet(dates.departure, i_dep);

[r_earth, v_earth] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
    kep_dep(4), kep_dep(5), kep_dep(6), mu_Sun);

time_dep = datetime( mjd20002date( dates.departure ) );

% asteroids parameters

[kep_NEO, ~] = ephNEO( dates.arrival, i_NEO);

[r_NEO, v_NEO] = kep2car(kep_NEO(1), kep_NEO(2), kep_NEO(3),...
    kep_NEO(4), kep_NEO(5), kep_NEO(6), mu_Sun);

time_arr = datetime( mjd20002date( dates.arrival ) );

% Vinf
TOF1 = seconds(diff([time_dep; time_fb]));
TOF2 = seconds(diff([time_fb; time_arr]));

[~,~,~,~,~,V_before_fb,~,~] = lambertMR(r_earth, r_saturn,TOF1, mu_Sun, 0, 0, 2);
[~,~,~,~,V_after_fb,~,~,~] = lambertMR(r_saturn, r_NEO,TOF2, mu_Sun, 0, 0, 2);

vinf_m = V_before_fb'-v_saturn;
vinf_p = V_after_fb'-v_saturn;

% normal to the manoeuvre plane
u = cross(vinf_m,vinf_p);
u = u/norm(u);

% Turning angle
delta = acos( vinf_m'*vinf_p / ( norm(vinf_m)*norm(vinf_p) ) );

% Radius of pericentre
fun = @(rp) 1e4* ( asin( 1 / (1+( rp*norm(vinf_m)^2) / mu_Sat) ) + ...
    asin(1/(1+(rp*norm(vinf_p)^2)/mu_Sat)) - delta );
R_Sat = astroConstants(26);
options = optimoptions('fsolve','OptimalityTolerance',1e-8);
rp=fsolve(fun, R_Sat,options);

if rp <= R_Sat || rp > R_Sat*906.9
    error('Unfeasible radius of pericentre for fly-by maneouvre')
end

% pericentre's coordinates
delta_m = asin( 1 / (1+( rp*norm(vinf_m)^2) / mu_Sat) );
rotation_angle = delta_m - pi/2;

v_rot_rp = vinf_m*cos(rotation_angle) + cross(u,vinf_m) * sin(rotation_angle) ...
    + u*dot(u,vinf_m)*(1-cos(rotation_angle));

rp = rp * v_rot_rp/norm(v_rot_rp);

v_rot_v = vinf_m*cos(rotation_angle + pi/2) + cross(u,vinf_m) * sin(rotation_angle+pi/2) ...
    + u*dot(u,vinf_m)*(1-cos(rotation_angle+pi/2));
vp_norm = v_rot_v/norm(v_rot_v);
% Powered deltaV

vp_minus = sqrt(norm(vinf_m)^2+2*mu_Sat/norm(rp));
vp_plus = sqrt(norm(vinf_p)^2+2*mu_Sat/norm(rp));

delta_vp = abs (vp_plus - vp_minus);


%% Hyperbolic trajectory plot

tof=1e9;
settings.perturbations = false;
SOI = R_Sat*906.9;
settings.mu = mu_Sat;
s1_i = [rp; vp_minus*vp_norm];
s1_f = [rp; vp_plus*vp_norm]; 

options = odeset('RelTol', 1e-10,'AbsTol',1e-11,'Events', @event_SOI );
[ti, Yi] = ode113(@pert_tbp,[tof 0], s1_i,options,settings,SOI);
[tf, Yf] = ode113(@pert_tbp,[0 tof], s1_f,options,settings,SOI);

% total flyby duration
flyby_duration = ( ti(end)+tf(end) ) / 3600 /24
figure()
hold on;
grid on;

plot3(Yi(:, 1), Yi(:, 2), Yi(:, 3), 'LineWidth',2);
plot3(Yf(:, 1), Yf(:, 2), Yf(:, 3), 'LineWidth',2);

opt.Units = 'km';
planet3D('Saturn',opt);




