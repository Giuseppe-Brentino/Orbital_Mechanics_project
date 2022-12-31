function dY = pert_tbp(~,s,settings, drag)
%
% PROTOTYPE:
%   dY = pert_tbp(~,s,settings, drag);
%
% DESCRIPTION:
%   ODE system to solve the two-body problem (Keplerian motion)
%
% INPUT:
%   t [1]       time (can be omitted)  [T]
%   s [6,1]     state of the body (position [L] and velocity [L/T]) [x y z vx vy vz]
%
%   settings    struct containing the following parameters: 
%                 - settings.perturbations  true or false, activates
%                                           perturbations
%                 - settings.mu [1]         planetary constant [km^3/s^2]
%
% (the following inpunts are required only if settings.perturbation = true)
%
%                 - settings.J2E     [1] gravitational field constant [-]
%                 - settings.RE      [1] planet's radius [km]
%                 - settings.w_E     [1] planet's angular velocity[rad/s]
% 
%   drag        struct containing the following parameters:
%                 - drag.Area_mass [1] reference area over mass [m^2/kg]
%                 - drag.c_d drag  [1] coefficient [-] 
%
% OUTPUT:
%   dY [6,1]    vector containing the derivative of the states (velocity[L/T] 
%               and acceleration [L/T^2]) [vx vy vz ax ay az]
%
% FUNCTIONS CALLED:
%   a_pert.m
%
% AUTHORS:
%   Giuseppe Brentino, Virginia Di Biagio Missaglia, Roberto Pistone
%   Nascone, 2022
%--------------------------------------------------------------------------

x = s(1);
y = s(2);
z = s(3);
vx = s(4);
vy = s(5);
vz = s(6);

v = [vx vy vz]';        % velocity vector
p = [x y z]';           % position vector

perturbations = settings.perturbations;

mu = settings.mu;

if perturbations
    settings.ref_sys = 'CAR';
    a_p = a_pert([p; v], settings, drag);
else
    a_p = zeros(3, 1);
end

r = norm(p);
dv = (-mu/(r^3) * p) + a_p;     % acceleration vector


dY(1:3) = v;
dY(4:6) = dv;

dY = dY';
