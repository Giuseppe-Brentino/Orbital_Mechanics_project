function dkep = GaussEq_rsw(~, kep, settings, drag)
%
% Gauss' equations in RSW (radial, transversal, out-of-plane) reference
% system for perturbed two-body problem
%
% INPUT:
%   t [1]       time (can be omitted)  [s]
%
%   s [6,1]     keplerian elements [a, e, i, OM, om, true anomaly] [km,rad]
%
%   settings    struct containing the following parameters: 
%                   - settings.mu      [1] planetary constant [km^3/s^2]
%                   - settings.J2E     [1] gravitational field constant [-]
%                   - settings.RE      [1] planet's radius [km]
%                   - settings.w_E     [1] planet's angular velocity[rad/s]
% 
%   drag        struct containing the following parameters:
%                   - drag.Area_mass [1] reference area over mass [m^2/kg]
%                   - drag.c_d drag  [1] coefficient [-] 
%
% OUTPUT:
%   dkep [6,1]  vector containing the derivative of the keplerian elements
%
% FUNCTIONS CALLED:
%   a_pert.m
%
% AUTHORS:
%   Giuseppe Brentino, Virginia Di Biagio Missaglia, Nicolò Galletta
%   Roberto Pistone Nascone, 2022
%--------------------------------------------------------------------------

a  = kep(1);
e  = kep(2);
i  = kep(3);
om = kep(5);
theta = kep(6);

mu = settings.mu;

p = a * (1 - e^2);
r = p / (1 + e*cos(theta));
h = sqrt(p * mu);

settings.ref_sys = 'RSW';
a_p_RSW = a_pert(kep, settings, drag);

a_r=a_p_RSW(1);
a_s=a_p_RSW(2);
a_w=a_p_RSW(3);


da = 2* a^2 / h * ( e * sin(theta)* a_r + p/r * a_s);

de = 1/h * ( p * sin(theta)* a_r+ ((p + r) * cos(theta)+ r * e) *a_s );

di = r * cos(theta+ om) / h * a_w;

dOM = r * sin(theta+ om) / h * a_w / sin(i);

dom = 1/h/e * ( -p * cos(theta) * a_r + (p + r) * sin(theta)* a_s)...
    - r * sin(theta + om) *cos(i ) / h / sin(i) * a_w;

dtheta = h/r^2 + 1/h/e * ( p * cos(theta) * a_r - (p + r) * sin(theta)* a_s);


dkep = [da, de, di, dOM, dom, dtheta]';

