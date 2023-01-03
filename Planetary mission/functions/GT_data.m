function [lon,lat] = GT_data(T,keplerian,settings,drag,theta_G0,t0)
%
% PROTOTYPE:
%  [lon,lat] = GT_data(T,t_span,keplerian,settings,drag,theta_G0,t0);
% 
% DESCRIPTION:
%   This function propagates the orbit for the given input and returns the 
%   longitude and the latitude for each timestep of the propagation.
%
% INPUT:
%   T [1]       Total period in which propagate the orbit   [s]
%
%   keplerian   struct containing the following parameters:
%                   - keplerian.a           [1] semi-major axis [km]
%                   - keplerian.e           [1] eccentricity [-]
%                   - keplerian.i           [1] inclination [rad]
%                   - keplerian.OM          [1] RAAN [rad]
%                   - keplerian.om          [1] argument of pericentre [rad]
%                   - keplerian.theta       [1] true anomaly [rad]
%
%   settings    struct containing the following parameters: 
%                   - settings.mu           [1] planetary constant [km^3/s^2]
%                   - settings.J2E          [1] gravitational field constant [-]
%                   - settings.RE           [1] planet's radius [km]
%                   - settings.w_E          [1] planet's angular velocity[rad/s]
%                   - settings.perturbations    true or false, activates perturbations
%                                               
%   drag        struct containing the following parameters:
%                   - drag.Area_mass [1] reference area over mass [m^2/kg]
%                   - drag.c_d drag  [1] coefficient [-] 
%
%  theta_G0[1]  Greenwich meridian's longitude @ t0 [rad]
%
%   t0[1]       Initial time instant to compute the groundtrack [s]
%
% OUTPUT:
%   lon [?,1]  vector containing the longitude for each timestep       [deg]
%   lat [?,1]  vector containig the latitude for each timestep         [deg]
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

% orbit propagation with to ODE113
[r0, v0] = kep2car(a, e, i, OM, om, theta, mu);
s0 = [r0; v0];                              

options = odeset('RelTol', 1e-13,'AbsTol',1e-14 );
[t, Y] = ode113(@pert_tbp, [0 T], s0, options, settings, drag);

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

