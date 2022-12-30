function [lon,lat] = GT_data(T,t_span,keplerian,settings,drag,theta_G0,t0)
%
% This function propagates the orbit for the given input and returns the 
% longitude and the latitude for each timestep of the propagation
%
% INPUT:
%   T [1]       time  [s]
%
%   t_span []     
%
%   keplerian
%
%   settings    struct containing the following parameters: 
%                   - settings.mu           [1] planetary constant [km^3/s^2]
%                   - settings.J2E          [1] gravitational field constant [-]
%                   - settings.RE           [1] planet's radius [km]
%                   - settings.w_E          [1] planet's angular velocity[rad/s]
%                   - settings.GT_ratio     [1] satellite's revs wrt planet's revs [-]
%                   - settings.perturbations    true or false, activates perturbations
%                                               
%   drag        struct containing the following parameters:
%                   - drag.Area_mass [1] reference area over mass [m^2/kg]
%                   - drag.c_d drag  [1] coefficient [-] 
%
%  theta_G0[1]
%
%   t0[1]
%
% OUTPUT:
%   lon []  vector containing the longitude for each timestep       [deg]
%   lat []  vector containig the latitude for each timestep         [deg]
%  
% FUNCTIONS CALLED:
%   kep2car.m, pert_tbp.m, groundTrack.m
%
% AUTHORS:
%   Giuseppe Brentino, Virginia Di Biagio Missaglia, Nicolò Galletta
%   Roberto Pistone Nascone, 2022
%--------------------------------------------------------------------------

a = keplerian.a;
e = keplerian.e;
i = keplerian.i;
OM = keplerian.OM;
om = keplerian.om;
theta = keplerian.theta;

mu = settings.mu;
w_E = settings.w_E;
GT_ratio = settings.GT_ratio;

orbit = 2*pi*sqrt(keplerian.a^3/settings.mu);
t_sample = 0:orbit/t_span:T;

% compute the new semi-major axis if the period is the modified one
if T == (2*pi/w_E)/GT_ratio
    a = (mu*(T/(2*pi))^2 )^(1/3);
    orbit = (2*pi/w_E)/GT_ratio;
    t_sample = 0:orbit/t_span:(23*3600 + 56*60 + 4)*4;
end

if settings.perturbations
    if T == (2*pi/w_E)/GT_ratio
        a = (mu*(T/(2*pi))^2 )^(1/3);
        orbit = (2*pi/w_E)/GT_ratio;
        t_sample = 0:orbit/t_span:(23*3600 + 56*60 + 4);
    end
end

% orbit propagation with to ODE113
[r0, v0] = kep2car(a, e, i, OM, om, theta, mu);
s0 = [r0; v0];                              % initial state vector

options = odeset('RelTol', 1e-13,'AbsTol',1e-14 );
[t, Y] = ode113(@pert_tbp, t_sample, s0, options, settings, drag);

% longitude and latitude computation
[lon, lat] = groundTrack (Y, theta_G0, t, w_E, t0);

for j = 2:length(lon)
    if lon(j-1)*lon(j) < 0 && abs(lon(j)) > 10
        lon = [lon(1:j-1); NaN; lon(j:end)];
        lat = [lat(1:j-1); NaN; lat(j:end)];
   end

    if lat(j-1)*lat(j) < 0 && abs(lat(j)) > 10
        lon = [lon(1:j-1); NaN; lon(j:end)];
        lat = [lat(1:j-1); NaN; lat(j:end)];
    end
end

