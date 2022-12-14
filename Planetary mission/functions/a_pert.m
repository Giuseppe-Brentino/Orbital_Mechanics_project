function a_p = a_pert(car_kep, settings, drag)
%
% PROTOTYPE:
%  a_p = a_pert(car_kep, settings, drag);
%
% DESCRIPTION:
%   The function computes the total perturbing accelerations due to planet
%   oblateness(J2) and to the aerodynamic drag of the Earth's atmosphere. 
%   For the aerodynamic drag in the function an exponential atmospheric 
%   model is implemented based on David A. Vallado - "Fundamentals of 
%   Astrodynamics and Applications" (1997). Atmospheric density model is 
%   valid only for planet Earth.
%
% INPUT :
%	car_kep [6]     state vector in cartesian [x y z vx vy vz] [km,km/s] 
%                   or keplerian elements [a e i Om om theta] [km, rad]
% 
%   settings        struct containing the following parameters: 
%                   - settings.mu      [1] planetary constant [km^3/s^2]
%                   - settings.J2E     [1] gravitational field constant [-]
%                   - settings.RE      [1] planet's radius [km]
%                   - settings.w_E     [1] planet's angular velocity[rad/s]
%                   - settings.ref_sys [1] equal to 'CAR' or 'RSW'
% 
%   drag            struct containing the following parameters:
%                   - drag.Area_mass [1] reference area over mass [m^2/kg]
%                   - drag.c_d drag  [1] coefficient [-] 
%
% OUTPUT:
%	a_p [3,1]    	Perturbing acceleration
%
%   Note: It is possible to compute the acceleration in Cartesian 
%   or RSW (radial-transversal) reference frame by specifying 
%   settings.ref_sys = 'CAR' or 'RSW'.
%
% FUNCTIONS CALLED:
%   (none)
%
% AUTHORS:
%   Virginia Di Biagio Missaglia, Nicolò Galletta, 
%   Roberto Piscone Nascone, 2022
% -------------------------------------------------------------------------

Area_mass = drag.Area_mass;
c_d = drag.c_d;

mu = settings.mu;
J2 = settings.J2E;
R_planet = settings.RE;
w_planet = settings.w_E;
ref_sys = settings.ref_sys;

switch ref_sys
    case 'CAR'
        x = car_kep(1);
        y = car_kep(2);
        z = car_kep(3);
        r_car(1:3, 1) = car_kep(1:3);
        v_car(1:3, 1) = car_kep(4:6);
        r = norm(car_kep(1:3));

    case 'RSW'
        a  = car_kep(1);
        e  = car_kep(2);
        i  = car_kep(3);
        OM = car_kep(4);
        om = car_kep(5);
        theta = car_kep(6);

        p = a * (1 - e^2);
        r = p / (1 + e*cos(theta));
        h = sqrt(p * mu);

        v_rsw = mu/h* [e*sin(theta); 1+e*cos(theta); 0];

        R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
        R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
        R3_omplusth = [cos(om+theta) sin(om+theta) 0;
                      -sin(om+theta) cos(om+theta) 0;
                      0              0             1];        
end

% Exponential atmospheric model
h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 ...
    300 350 400 450 500 600 700 800 900 1000]';                                 % initial base altitude [km]

h0_lim = 1200;                                                                  % maximum altitude at which aerodynamic drag has effect [km]

rho0_vect = [1.225, 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3, 3.206*1e-4 ...
    8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 2.438*1e-8 ...
    8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 2.789*1e-10 7.248*1e-11 ...
    2.418*1e-11 9.158*1e-12 3.725*1e-12 1.585*1e-12 6.967*1e-13 1.454*1e-13...
    3.614*1e-14 1.170*1e-14 5.245*1e-15 3.019*1e-15]' * 1e3;                    % initial nominal density [kg/m^2/Km]

H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
    7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 53.298 ...
    58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';                  % initial scale height [km]

h_alt = r - R_planet;

if h_alt <= 0
    error('Altitude lower than 0 km');
end

if h_alt >= h0_lim
    rho = 0;
elseif h_alt >= h0_vect(end)
    h0 = h0_vect(end);
    H = H_vect(end);
    rho0 = rho0_vect(end);
    rho = rho0 * exp(-(h_alt-h0)/H);
else
    for j = 1:(length(h0_vect)-1)
        if h_alt >= h0_vect(j) && h_alt < h0_vect(j+1)
            rho0 = rho0_vect(j);
            H = H_vect(j);
            h0 = h0_vect(j);
            rho = rho0 * exp(-(h_alt-h0)/H);
        end
    end
end

% Perturbing accelerations
switch ref_sys
    case 'CAR'
        % J2 perturbation
        c = 3/2 * (J2*mu*R_planet^2)/r^4;
        a_J2 = c * [ x/r * ( 5*(z/r)^2 - 1 ); ...
            y/r * ( 5*(z/r)^2 - 1 ); ...
            z/r * ( 5*(z/r)^2 - 3 ) ] ;

        % Drag perturbation
        w_E_car = w_planet * [0 0 1]';
        v_rel = v_car - cross(w_E_car, r_car);
        v_rel_norm = norm(v_rel);

        a_drag = -0.5 * (Area_mass*c_d) * rho * v_rel_norm^2 * (v_rel/v_rel_norm);

    case 'RSW'
        % J2 perturbation
        c = - 3/2 * (J2*mu*R_planet^2)/r^4;

        a_J2 = c * [ 1 - 3 * (sin(i))^2 * (sin(theta + om))^2; ...
            (sin(i))^2 * sin(2*(theta + om)); ...
            sin(2*i) * sin(theta + om)];
        
        % Drag perturbation
        w_E_rsw = R3_omplusth * R1_i * R3_OM * [0; 0; w_planet];
        v_rel_rsw = v_rsw - cross(w_E_rsw, [1;0;0]*r);
        v_rel_norm_rsw = norm(v_rel_rsw);

        a_drag = -0.5 * (Area_mass*c_d) * rho * v_rel_norm_rsw^2 * (v_rel_rsw/v_rel_norm_rsw);
        
end

a_p = a_drag + a_J2;
