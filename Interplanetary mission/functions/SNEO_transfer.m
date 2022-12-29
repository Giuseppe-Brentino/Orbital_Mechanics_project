function [delta_v_arr,delta_v_dep,VI,TPAR] = SNEO_transfer(t_dep, t_arr, mu, i_dep, i_arr)
%
% PROTOTYPE:
%  [delta_v_arr,delta_v_dep,VI,TPAR] = SNEO_transfer(t_dep, t_arr, mu,
%  i_dep, i_arr);
%
% DESCRIPTION:
%   The function computes the delta velocity needed: to perform the
%   injection manoeuvre (from initial orbit around Saturn to the Transfer
%   arc) and the arrival manoeuvre (from the Transfer arc to the final
%   orbit around the designed asteroid).
%
%  INPUT :
%	t_dep[1]    Departure time                              [MJD2000]
%	t_arr[1]    Arrival time                                [MJD2000]
%   mu[1]       Gravitational constant of the Sun           [km^3/s^2]
%   i_dep[1]    Index of the departing planet to be used in the function
%               uplanet
%   i_arr[1]    Index of the arriving asteroid to be used in the function
%               ephNEO
%
%  OUTPUT:
%	delta_v_arr[1]  Delta velocity needed for the arrival manoeuvre     [km/s]
%   delta_v_dep[1]  Delta velocity needed for the injection manoeuvre   [km/s]
%   VI[1x3]         Initial velocity vector from lambertMR          [km/s]
%   TPAR[1]         Parabolic ToF from lambertMR                    [s]
%
%  FUNCTIONS CALLED:
%   mjd20002date.m, uplanet.m, kep2car.m, ephNEO.m, lambertMR.m 
%
% AUTHORS:
%   Virginia Di Biagio Missaglia, Roberto Pistone Nascone, Giuseppe
%   Brentino, Nicolò Galletta
%
% -------------------------------------------------------------------------

%% Compute Time Of Fly of the transfer orbit
time_arr = datetime(mjd20002date(t_arr));
time_dep = datetime(mjd20002date(t_dep));

TOF = seconds(diff([time_dep; time_arr]));
if TOF >= 0

    [kep_dep, ~] = uplanet(t_dep, i_dep);

    [r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
        kep_dep(4), kep_dep(5), kep_dep(6), mu);

    [kep_arr, ~] = ephNEO(t_arr, i_arr);

    [r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
        kep_arr(4), kep_arr(5), kep_arr(6), mu);

    [~,~,~,~,VI,VF,TPAR,~] = lambertMR(r_dep, r_arr, TOF, mu, 0, 0, 2);
    delta_v_dep = norm(VI - v_dep');
    delta_v_arr = norm(VF - v_arr');

else
    delta_v_dep = NaN;
    delta_v_arr = NaN;
    TPAR = NaN;
    VI = NaN;
end

