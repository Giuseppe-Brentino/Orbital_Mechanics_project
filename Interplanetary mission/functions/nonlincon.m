function [C,Ceq] = nonlincon(v_inf, t_arr, t_dep, i_dep, i_arr, mu)
Ceq = [];

time_arr = datetime(mjd20002date(t_arr));
time_dep = datetime(mjd20002date(t_dep));

TOF = seconds(diff([time_dep; time_arr]));

[kep_dep, ~] = uplanet(t_dep, i_dep);

[r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
    kep_dep(4), kep_dep(5), kep_dep(6), mu);

[kep_arr, ~] = uplanet(t_arr, i_arr);

[r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
    kep_arr(4), kep_arr(5), kep_arr(6), mu);


[~,~,~,~,VI,~,~,~] = lambertMR(r_dep, r_arr,TOF, mu, 0, 0, 2);
C = norm(VI - v_dep') - v_inf;

end

