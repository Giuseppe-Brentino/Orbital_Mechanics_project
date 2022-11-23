function [departure_times, arrival_times, r_dep, v_dep, r_arr, v_arr] = timeWindowMj2000...
    (departure_date_e, departure_date_l, arrival_date_e, arrival_date_l, i_dep, i_arr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




t_dep_e = datetime(departure_date_e);
t_dep_l = datetime(departure_date_l);


t_arr_e = datetime(arrival_date_e);
t_arr_l = datetime(arrival_date_l);


t_dep = t_dep_e : days (1) : t_dep_l;
t_arr = t_arr_e : days (1) : t_arr_l;

dep_times = datevec(t_dep);
arr_times = datevec(t_arr);

dep_times = [dep_times; departure_date_l] ;

arr_times = [arr_times; arrival_date_l];

departure_times = zeros (size(dep_times, 1), 1);
arrival_times = zeros (size(arr_times, 1), 1);
r_dep = zeros(size(dep_times,1), 3);
r_arr = zeros(size(arr_times,1), 3);

v_dep = zeros(size(dep_times,1), 3);
v_arr = zeros(size(arr_times,1), 3);


for i = 1 : length(departure_times)

    departure_times(i) = date2mjd2000(dep_times(i, :));

    [kep_dep,mu_sun] = uplanet(departure_times(i), i_dep);

    [r_dep(i, :), v_dep(i, :)] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
        kep_dep(4), kep_dep(5), kep_dep(6), mu_sun);

end

for i = 1 : length(arrival_times)

    arrival_times(i) = date2mjd2000(arr_times(i, :));

    [kep_arr,mu_sun] = uplanet(arrival_times(i), i_arr);

    [r_arr(i, :),v_arr(i, :)] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
        kep_arr(4), kep_arr(5), kep_arr(6), mu_sun);

end








