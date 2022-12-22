function a_p = a_pert(settings, p, v, drag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
%
% INPUT:
%   
%
% OUTPUT:
%   
%
% Contributors:
%   Virginia Di Biagio Missaglia, Roberto Pistone Nascone, Nicolò Galletta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Area_mass = drag.Area_mass;
c_d = drag.c_d;

mu = settings.mu;
J2 = settings.J2E;
RE = settings.RE;
w_E = settings.w_E;

x = p(1);
y = p(2);
z = p(3);
vx = v(1);
vy = v(2);
vz = v(3);


r = norm(p);


h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 ...
    300 350 400 450 500 600 700 800 900 1000]';                                             % initial base altitude [km]

rho0_vect = [1.225, 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3, 3.206*1e-4 ...
    8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 2.438*1e-8 ...
    8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 2.789*1e-10 7.248*1e-11 ...
    2.418*1e-11 9.158*1e-12 3.725*1e-12 1.585*1e-12 6.967*1e-13 1.454*1e-13...
    3.614*1e-14 1.170*1e-14 5.245*1e-15 3.019*1e-15]';                                      % initial nominal density [kg/m^3]

H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
    7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 53.298 ...
    58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';                              % initial scale height [km]

w_E_vect = w_E * [0 0 1]';
v_rel = v - cross(w_E_vect, p);          % initial relative velocity vector [km/s]
v_rel_norm = norm(v_rel);               % initial relative velocity vector's norm [km/s]



% J2 perturbation
c = 3/2 * (J2*mu*RE^2)/r^4;
a_J2 = c * [ x/r * ( 5*(z/r)^2 - 1 ); ...
    y/r * ( 5*(z/r)^2 - 1 ); ...
    z/r * ( 5*(z/r)^2 - 3 ) ] ;


% Drag perturbation
h = r - RE;

if h < 0
    error('Altitude lower than 0 km');
end

for j = 1:length(h0_vect)
    if h >= 1000
        h0 = h0_vect(end);
        H = H_vect(end);
        rho0 = rho0_vect(end);

    elseif h >= h0_vect(j) && h < h0_vect(j+1)
        rho0 = rho0_vect(j);
        H = H_vect(j);
        h0 = h0_vect(j);

    end
end


rho = rho0 * exp(-(h-h0)/H);

a_drag = -0.5 * (Area_mass*c_d) * rho * v_rel_norm^2 * (v_rel/v_rel_norm);

% perturbated acceletation
a_p = a_J2 + a_drag;
