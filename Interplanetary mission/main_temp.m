clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');


mu_Sun = astroConstants(4);
mu_E = astroConstants(13);
mu_Sat = astroConstants(16);




%% Saturn-NEO Leg

departure_date_e = [2031, 8, 2, 0, 0, 0]; 
departure_date_l = [2064, 1, 28, 23, 59, 59];

arrival_date_e = [2032, 8, 2, 0, 0, 0];
arrival_date_l = [2065, 1, 28, 23, 59, 59];

delta_t = days(90);

[t_fb, t_arr] = timeWindowMj2000...
    (departure_date_e, departure_date_l, arrival_date_e, arrival_date_l, delta_t);

i_dep_1 = 6;  % Saturn
i_arr_1 = 21; % NEO


m = length(t_fb);
n = length(t_arr);

% indeces i, m for t_fb_1 elements
% indeces j, n for t_arr elements

delta_v_arr = zeros(m, n);

VF_1 = cell(m, n);
VI_1 = cell(m, n);

T_par_1 = zeros(m, n);
for i = 1:m
    for j = 1:n
        [delta_v_arr(i, j), VF_1{i, j}, VI_1{i,j}, T_par_1(i, j)] = SNEO_transfer (t_fb(i), t_arr(j), mu_Sun, i_dep_1, i_arr_1);
      
    end
end


%% Cheapest departure/arrival to NEO dates 

t_fb_1_logic = false(length(t_fb), length(t_arr));

i = 1;
 for k = 1:size(delta_v_arr, 1)
    for j = 1:size(delta_v_arr, 2)
        time_fb_1 = datetime(mjd20002date(t_arr(j)));
         time_dep_1 = datetime(mjd20002date(t_fb(k)));

        if delta_v_arr(k, j) < 3.5 && T_par_1(k, j) < ( seconds(diff([time_dep_1; time_fb_1])) )
            
            t_fb_1_logic(k, j) = true;

            S_N_data{i} = [t_fb(k); t_arr(j); delta_v_arr(k,j);VI_1{k,j}';VF_1{k,j}];
            i = i+1;

        end
    end
 end


%% PorkChop plot Saturn-NEO

figure()
hold on;
contour(t_fb, t_arr, delta_v_arr', [0:0.5:8], 'LineWidth',2);
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([0, 8]);

%% Comparison
fb_best_times = find(~all(t_fb_1_logic==0,2)); % find the rows with min delta v (fb date)

t_fb = t_fb(fb_best_times);

delta_v_arr = delta_v_arr(fb_best_times,:); % extract delta_v rows of best fb date

%% Earth-Saturn Leg

i_dep = 3; % Earth
i_arr = 6; % Saturn

t_dep = linspace(t_fb(1)-730, t_fb(end)-180, 100);

m = length(t_dep);
n = length(t_fb);

% indeces i, m for t_dep elements
% indeces j, n for t_fb elements

delta_v_dep = zeros(m, n);
VF = cell(m, n);
VI = cell(m, n);
T_par = zeros(m, n);
for i = 1:m
    for j = 1:n
        [delta_v_dep(i, j), VF{i, j}, VI{i, j}, T_par(i, j)] = ES_transfer (t_dep(i), t_fb(j), mu_Sun, i_dep, i_arr);
      
    end
end


%% PorkChop plot Saturn-Earth

figure()
hold on;
contour(t_dep, t_fb, delta_v_dep', 10:2:30, 'LineWidth',2);
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([5, 30]);

%% Cheapest departure/arrival to NEO dates 

i = 1;
for k = 1:size(delta_v_dep, 1)
    for j = 1:size(delta_v_dep, 2)
 
        time_fb = datetime(mjd20002date(t_fb(j)));
        time_dep = datetime(mjd20002date(t_dep(k)));

        if delta_v_dep(k, j) < 25 && T_par(k, j) < ( seconds(diff([time_dep; time_fb])) )
            E_S_data{i} = [t_dep(k); t_fb(j); delta_v_dep(k,j); VI{k,j}; VF{k,j}'];
            i = i+1;

        end
    end
end

%% Find best date
k = 1;
for i = 1 : length(S_N_data)
    for j = 1 : length(E_S_data)
        if S_N_data{i}(1) == E_S_data{j}(2)
            delta_v = E_S_data{j}(4) + norm(E_S_data{j}(5:7)-S_N_data{i}(4:6))...
                + S_N_data{i}(7);
            mission_data(:,k) =[E_S_data{j}(1); E_S_data{j}(2); S_N_data{i}(2); delta_v];
            k = k+1;
        end
    end
end

[balue,index] = min(mission_data(4,:));

%% Optimization

guess = [mission_data(1,index);mission_data(2,index);mission_data(3,index)];

lb = guess - 180; 
ub = guess + 180; 
 
options = optimoptions("fmincon", "Display","iter",'Algorithm','interior-point');
[u, DeltaV_opt, EXITFLAG, OUTPUT] = fmincon(@(u)obj_Fcn(u,i_dep,i_arr,i_arr_1),...
    guess, [], [], [], [], lb, ub,[], options);

%% Saturn Fly-by










