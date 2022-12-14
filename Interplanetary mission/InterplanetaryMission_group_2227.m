clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');


mu_Sun = astroConstants(4);
mu_E = astroConstants(13);
mu_Sat = astroConstants(16);

%% Earth-Saturn Leg

% in italiano: stiamo considerando che metà del tempo della missione viene
% impiegato per arrivare a S mentre la restante parte per arrivare a Neo
% PROVARE CON FINESTRA COMPLETA ANCHE T-S 
departure_date_e = [2030, 8, 2, 0, 0, 0];
departure_date_l = [2040, 1, 28, 23, 59, 59]; 


arrival_date_e = [2031, 8, 2, 0, 0, 0]; 
arrival_date_l = [2050, 1, 28, 23, 59, 59]; 


i_dep = 3; % Earth
i_arr = 6; % Saturn

[t_dep, t_fb] = timeWindowMj2000...
    (departure_date_e, departure_date_l, arrival_date_e, arrival_date_l);


m = length(t_dep);
n = length(t_fb);

% indeces i, m for t_dep elements
% indeces j, n for t_arr elements

delta_v_dep = zeros(m, n);
VF = cell(m, n);
T_par = zeros(m, n);
for i = 1:m
    for j = 1:n
        [delta_v_dep(i, j), VF{i, j}, T_par(i, j)] = ES_transfer (t_dep(i), t_fb(j), mu_Sun, i_dep, i_arr);
      
    end
end



%% PorkChop plot

figure()
hold on;
contour(t_dep, t_fb, delta_v_dep', [10:2:30], 'LineWidth',2);
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([5, 30]);

%% Cheapest departure/arrival to Sat dates 

t_dep_logic = false(length(t_dep), length(t_fb));
t_fb_logic = false(1, length(t_fb));

for k = 1:size(delta_v_dep, 1)
    for j = 1:size(delta_v_dep, 2)
        time_fb = datetime(mjd20002date(t_fb(j)));
        time_dep = datetime(mjd20002date(t_dep(k)));

        if delta_v_dep(k, j) < 20 && T_par(k, j) < ( seconds(diff([time_dep; time_fb])) )
            t_fb_logic(j) = true;
            [~, index] = min(delta_v_dep(:, j));
            t_dep_logic(index, j) = true;

        end
    end
end


