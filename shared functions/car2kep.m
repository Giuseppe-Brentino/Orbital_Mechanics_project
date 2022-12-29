function [a, e, i, OM, om, theta] = car2kep(r, v, mu)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Conversion from Keplerian elements to Cartesian coordinates. Angles in
%  radians
%
% OUTPUT:
%   a   [1x1] Semi-major axis           [km]
%   e   [1x1] Eccentricity              [-]
%   i   [1x1] Inclination               [rad]
%   OM  [1x1] RAAN                      [rad]
%   om  [1x1] Pericentre anomaly        [rad]
%   th  [1x1] True anomaly              [rad]
%   mu  [1x1] Gravitational parameter   [km^3/s^2]
%
% INPUT:
%   r   [3x1] Position vector           [km]
%   v   [3x1] Velocity vector           [km/s]
%
% Contributors:
%   Nicolò Galletta, Virginia di Biagio Missaglia, Giuseppe Brentino,
%   Roberto Pistone Nascone
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                                            % reference system unit vectors
I = [1 0 0]';
J = [0 1 0]';
K = [0 0 1]';


r_norm = norm(r); % magnitude of position vector[km]
 
v_norm = norm(v);                                                           % magnitude of veocity vector [km/s]

a = 1 / (2 / r_norm - v_norm^2 / mu);                                       % semi-major axis [km]

h_vect = cross(r,v);                                                        % specific angular momentum vector [km^2/s]
h = norm(h_vect);

e_vect = 1 / mu * ( cross(v,h_vect) - mu * r / r_norm );                    % eccentricity vector [-]
e = norm(e_vect);
 
i = acos ( dot(h_vect,K) / h );    % orbit's inclination [rad]

%Check on the inclination

if i==0
  N_vect = I;
  n_vect = I; 
else
N_vect = cross(K,h_vect);
n_vect = cross(K,h_vect)/norm(cross(K,h_vect)); % Node line
end 

%check on the eccentricity

if e < 1e-10
    e=0; 
    e_vect = N_vect; 
    om = 0; 
else
    om = acos( dot(n_vect,e_vect) / e );                                        % argument of pericentre [rad]

    if ( dot(e_vect,K) < 0 )

    om = 2*pi - om;

    end

end

%RAAN 

OM = acos ( dot(n_vect,I) ); 


if ( dot(n_vect,J) < 0 )
    OM = 2*pi - OM;
end


check = abs(dot(r,e_vect) - ( r_norm * e )) ; 
    
if check < 1e-8
    theta = 0; 
else 

    if e ==0

    theta = acos( dot(r,e_vect) / (r_norm) ); 

    else

    theta = acos( dot(r,e_vect) / ( r_norm * e ) ); 

    end 

end 

if (dot(v,r) < -1e-8)
    theta = 2*pi - theta;
end
  
