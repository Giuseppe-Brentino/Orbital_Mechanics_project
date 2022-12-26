function [delta_v, VF,delta_v_dep, TPAR] = ES_transfer(t_dep, t_arr, mu, i_dep, i_arr)
%{
INPUT:
  t_dep [1] departure time                                     [MJD2000]
  t_arr [1] arrival time                                       [MJD2000]
  mu    [1] gravitational constant of the sun
  i_dep [1] index of the departing planet to be used in the function
            uplanet
  i_arr [1] index of the arriving planet to be used in the function
            uplanet
OUTPUT:
  delta_v [?, 1] vectors of norm of every possible delta velocity [km/s]

Contributors:
  Nicolò galletta, Roberto Pistone Nascone, Giuseppe
  Brentino

%} 


time_arr = datetime(mjd20002date(t_arr));
time_dep = datetime(mjd20002date(t_dep));

TOF = seconds(diff([time_dep; time_arr]));

if TOF >= 0

    [kep_dep, ~] = uplanet(t_dep, i_dep);

    [r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
        kep_dep(4), kep_dep(5), kep_dep(6), mu);

    [kep_arr, ~] = uplanet(t_arr, i_arr);

    [r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
        kep_arr(4), kep_arr(5), kep_arr(6), mu);


    [~,~,~,~,VI,VF,TPAR,~] = lambertMR(r_dep, r_arr,TOF, mu, 0, 0, 2);
    
    delta_v = norm(VI - v_dep') + norm(VF - v_arr');
    delta_v_dep = norm(VI - v_dep');

else
    delta_v = NaN;
    VF = NaN;
    TPAR = NaN;
    delta_v_dep = NaN;
end


