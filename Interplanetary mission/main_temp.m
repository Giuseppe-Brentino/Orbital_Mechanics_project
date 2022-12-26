%% Orbital Mechanics Project - Interplanetary mission

clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
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

        if S_N.deltaV_arr(k, j) < 25 && S_N.Tpar(k, j) < ( seconds(diff([temp.time_dep; temp.time_arr] ) ) ) % S_N.deltaV prima era S_N.deltaV_arr
            
            dates.flyby_logic(k, j) = true;

            S_N.data{i} = [dates.flyby(k); dates.arrival(j);...
                S_N.deltaV_arr(k,j); S_N.VI{k,j}'; S_N.deltaV_dep{k,j} ];

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
contour(dates.departure, dates.flyby, E_S.deltaV_dep', 10:2:30, 'LineWidth',2);
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([5, 30]);

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
 
options = optimoptions("fmincon", "Display","iter",'Algorithm','interior-point');
[u, DeltaV_opt, EXITFLAG, OUTPUT] = fmincon(@(u)obj_Fcn(u,i_dep,i_flyby,i_NEO),...
    guess, [], [], [], [], lb, ub,[], options);

dep_date = u(1);
fb_date  = u(2);
arr_date = u(3);

%% Saturn Fly-by

% saturn parameters
[kep_saturn,~] = uplanet( fb_date, i_arr);

[r_saturn, v_saturn] = kep2car(kep_saturn(1), kep_saturn(2), kep_saturn(3),...
    kep_saturn(4), kep_saturn(5), kep_saturn(6), mu_Sun);

time_fb = datetime( mjd20002date( fb_date ) );

% earth parameters
[kep_dep, ~] = uplanet(dep_date, i_dep);

[r_earth, v_earth] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
    kep_dep(4), kep_dep(5), kep_dep(6), mu_Sun);

time_dep = datetime( mjd20002date( dep_date ) );

% asteroids parameters

[kep_NEO, ~] = ephNEO( arr_date, i_arr_1);

[r_NEO, v_NEO] = kep2car(kep_NEO(1), kep_NEO(2), kep_NEO(3),...
    kep_NEO(4), kep_NEO(5), kep_NEO(6), mu_Sun);

time_arr = datetime( mjd20002date( arr_date ) );

% Vinf
TOF1 = seconds(diff([time_dep; time_fb]));
TOF2 = seconds(diff([time_fb; time_arr]));

[~,~,~,~,~,V_before_fb,~,~] = lambertMR(r_earth, r_saturn,TOF1, mu_Sun, 0, 0, 2);
[~,~,~,~,V_after_fb,~,~,~] = lambertMR(r_saturn, r_NEO,TOF2, mu_Sun, 0, 0, 2);

vinf_m = V_before_fb'-v_saturn;
vinf_p = V_after_fb'-v_saturn;

% Turning angle
delta = acos( vinf_m'*vinf_p / ( norm(vinf_m)*norm(vinf_p) ) );

% Radius of pericentre
fun = @(rp) asin( 1 / (1+( rp*norm(vinf_m)^2) / mu_Sat) ) + ...
    asin(1/(1+(rp*norm(vinf_p)^2)/mu_Sat)) - delta;
R_Sat = astroConstants(26);

rp=fsolve(fun, R_Sat);

if rp <= R_Sat || rp > R_Sat*906.9
    error('Unfeasible radius of pericentre for fly-by maneouvre')
end

% Powered deltaV

vp_minus = sqrt(norm(vinf_m)^2+2*mu_Sat/rp);
vp_plus = sqrt(norm(vinf_p)^2+2*mu_Sat/rp);

delta_vp = abs (vp_plus - vp_minus);

%% Hyperbolic trajectory plot

tof=1e8;
settings.perturbations = false;
SOI = R_Sat*906.9;
settings.mu = mu_Sat;
s1_i = [rp*[1 0 0]'; vp_minus*[0 1 0]'];
s1_f = [rp*[1 0 0]'; vp_plus*[0 1 0]']; 

options = odeset('RelTol', 1e-10,'AbsTol',1e-11,'Events', @event_SOI );
[ti, Yi] = ode113(@pert_tbp,[tof 0], s1_i,options,settings,SOI);
[tf, Yf] = ode113(@pert_tbp,[0 tof], s1_f,options,settings,SOI);

figure()
hold on;
grid on;
axis equal;
plot3(Yi(:, 1), Yi(:, 2), Yi(:, 3), 'LineWidth',2);
plot3(Yf(:, 1), Yf(:, 2), Yf(:, 3), 'LineWidth',2);
opt.Units = 'km';
planet3D('Saturn',opt);




