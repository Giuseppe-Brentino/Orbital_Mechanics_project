function [delta_v] = obj_Fcn(u,i_dep,i_fb,i_arr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ??
%
% INPUT:
%   u [3x1] departure, flyby and arrival dates  [MJD2000]
%
% OUTPUT:
%   delta_v [1] total delta v for the mission [km/s]
%
% Contributors:
%   Nicolò Galletta, Roberto Pistone Nascone, Giuseppe Brentino
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Earth data
[kep_dep, mu] = uplanet(u(1), i_dep);

[r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
    kep_dep(4), kep_dep(5), kep_dep(6), mu);

%% Saturn Data
[kep_fb, ~] = uplanet( u(2), i_fb);

[r_fb, v_fb] = kep2car(kep_fb(1), kep_fb(2), kep_fb(3),...
    kep_fb(4), kep_fb(5), kep_fb(6), mu);

%% First heliocentric leg parameters

time_dep = datetime( mjd20002date( u(1) ) );
time_fb = datetime( mjd20002date( u(2) ) );

TOF = seconds(diff([time_dep; time_fb]));

[~,~,~,~,VI,V_before_fb,~,~] = lambertMR(r_dep, r_fb,TOF, mu, 0, 0, 2);
delta_v1 = norm(VI - v_dep');

%%  Asteroids parameter

[kep_arr, ~] = ephNEO( u(3), i_arr);

[r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
    kep_arr(4), kep_arr(5), kep_arr(6), mu);

%% Second heliocentric leg parameters

time_arr = datetime( mjd20002date( u(3) ) );

TOF = seconds(diff([time_fb; time_arr]));

[~,~,~,~,V_after_fb,VF,~,~] = lambertMR(r_fb, r_arr,TOF, mu, 0, 0, 2);
delta_v2 = norm(v_arr' - VF);

%% DeltaV flyby

mu_Sat = astroConstants(16);

vinf_m = V_before_fb'-v_fb;
vinf_p = V_after_fb'-v_fb;

% Turning angle
delta = acos( vinf_m'*vinf_p / ( norm(vinf_m)*norm(vinf_p) ) );

% Radius of pericentre
fun = @(rp) 1e4* ( asin( 1 / (1+( rp*norm(vinf_m)^2) / mu_Sat) ) + ...
    asin(1/(1+(rp*norm(vinf_p)^2)/mu_Sat)) - delta );
R_Sat = astroConstants(26);


opt = optimset('Display','none', 'TolFun',1e-8);
rp = fsolve(fun, R_Sat,opt);

% Powered deltaV

vp_minus = sqrt(norm(vinf_m)^2+2*mu_Sat/rp);
vp_plus = sqrt(norm(vinf_p)^2+2*mu_Sat/rp);

delta_vp = abs (vp_plus - vp_minus);

%% Cost function

delta_v = delta_v1 + delta_v2 + delta_vp;
