function dY = pert_tbp (~,s,settings,varargin)
% ODE system to solve the two-body problem (Keplerian motion)
%
% INPUT:
%   t [1]: time (can be omitted)  [T]
%   s [6x1]: States of the body (position [L] and velocity [L/T]) [x y z vx vy vz]
%   settings: Structure containing the required simulation parameters, such
%             us earth radius and gravitational constant
%   varargin: [1x?] Cell that contains the name of the required
%   perturbations:
%   -  oblateness: takes into account the non-sphericity of the main body
%
% OUTPUT:
%   dY [6x1]: Vector containing the derivative of the states (velocity[L/T] 
%             and acceleration [L/T^2]) [vx vy vz ax ay az]
%
% Contributors:
%   Giuseppe Brentino, Virginia Di Biagio Missaglia, Roberto Pistone
%   Nascone
%

x = s(1);
y = s(2);
z = s(3);
vx = s(4);
vy = s(5);
vz = s(6);

v = [vx vy vz]';        % velocity vector
p = [x y z]';           % position vector

mu = settings.mu;

r = norm(p);            % distance between the two bodies
dv = -mu/(r^3) * p;     %accelration vector

if nargin == 4
    perturbations = varargin{1};

    for i = 1:length(perturbations)

        switch(perturbations{i})
            case 'oblateness'
                J2 = settings.J2E;
                R = settings.RE;
            
                c = 3/2 * (J2*mu*R^2)/r^4;
                a_J2 = c * [ x/r * ( 5*(z/r)^2 - 1 ); ...
                             y/r * ( 5*(z/r)^2 - 1 ); ...
                             z/r * ( 5*(z/r)^2 - 3 ) ] ;
                
                dv = dv + a_J2;
        end

    end

end
dY(1:3) = v;
dY(4:6) = dv;

dY = dY';
