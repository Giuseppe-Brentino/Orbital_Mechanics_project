clc; clearvars; close all;
 
addpath('../given functions/');
addpath('../given functions/time/');
addpath('../shared functions/');
addpath('./functions/');


mu_Sun = astroConstants(4);
mu_E = astroConstants(13);
mu_Sat = astroConstants(16);

v_inf = 3.5;
departure_date_e = [ 2030, 8, 2, 0, 0, 0]; 
departure_date_l = [ 2052, 1, 28, 23, 59, 59]; 


arrival_date_e = [2043, 8, 2, 0, 0, 0]; 
arrival_date_l = [2065, 1, 28, 23, 59, 59]; 


i_dep = 3; % Earth
i_arr = 6; % Saturn

[t_dep, t_arr] = timeWindowMj2000...
    (departure_date_e, departure_date_l, arrival_date_e, arrival_date_l);


m = length(t_dep);
n = length(t_arr);

% indeces i, m for t_dep elements
% indeces j, n for t_arr elements

delta_v = zeros(m, n);
for i = 1:m
    for j = 1:n
        delta_v(i, j) = obj_Fcn (t_dep(i), t_arr(j), mu_Sun, i_dep, i_arr);
      
    end
end

%% Hohmann ToF
[kep_E, ksun_E] = uplanet(date2mjd2000(departure_date_e), i_dep);
[kep_Sat, ksun_Sat] = uplanet(date2mjd2000(departure_date_e), i_arr);

r_E = kep_E(1);
r_Sat = kep_Sat(1);

a_H = (r_E + r_Sat)/2;
ToF_H = 2*pi*sqrt(a_H^3/ksun_E)/3600/24/365;


%% PorkChop Plot
% 
% for k = 1:size(delta_v, 1)
%     for j = 1:size(delta_v, 2)
%         if delta_v(k, j) > 20
%             delta_v(k, j) = NaN;
%         end
%     end
% end

figure()
hold on;
contour(t_dep, t_arr, delta_v', [10:5:30], 'LineWidth',2);
xlabel('Departure Time');
ylabel('Arrival Time');
colorbar();
clim([0, 30]);
