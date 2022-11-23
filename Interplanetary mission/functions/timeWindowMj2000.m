function [departure_times, arrival_times] = timeWindowMj2000...
    (departure_date_e, departure_date_l, arrival_date_e, arrival_date_l)

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

for i = 1 : length(departure_times)

    departure_times(i) = date2mjd2000(dep_times(i, :));

end

for i = 1 : length(arrival_times)

    arrival_times(i) = date2mjd2000(arr_times(i, :));

end










