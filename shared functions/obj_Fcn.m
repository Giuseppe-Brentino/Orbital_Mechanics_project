function [delta_v] = obj_Fcn(t_dep, t_arr, r_dep, r_arr, v1, v2, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ??
%
% INPUT:
%   t_dep [1] departure time                                     [MJD2000]
%   t_arr [1] arrival time                                       [MJD2000]
%   r_dep [?, 3] position vectors corresponding at departure time   [km]
%   r_arr [?, 3] position vectors corresponding at arrival time     [km]
%   v1    [?, 3] matrix velocity planet 1 for each departure time  [km/s]
%   v2 
%   mu
%
% OUTPUT:
%   delta_v [?, 1] vectors of norm of every possible delta velocity [km/s]
%
% Contributors:
%   Nicolò Galletta, Roberto Pistone Nascone
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_arr = mjd20002date(t_arr);
t_dep = mjd20002date(t_dep);

TOF = seconds(diff(t_dep, t_arr));

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r_dep, r_arr,TOF, mu, 0, 0, 2);
delta_v = norm(VI - v1) + norm(v2 - VF);


