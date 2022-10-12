function [r, v] = kep2car(a, e, i, OM, om, th, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v]= kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
%  Conversion from Keplerian elements to Cartesian coordinates. Angles in
%  radians
%  
% INPUT:
%   a   [1x1] Semi-major axis           [km]
%   e   [1x1] Eccentricity              [-]
%   i   [1x1] Inclination               [rad]
%   OM  [1x1] RAAN                      [rad]
%   om  [1x1] Pericentre anomaly        [rad]
%   th  [1x1] True anomaly              [rad]
%   mu  [1x1] Gravitational parameter   [km^3/s^2]
%
% OUTPUT:
%   r   [3x1] Position vector           [km]
%   v   [3x1] Velocity vector           [km/s]

if nargin==6
    mu=398600.433;
end

p = a*(1-e^2);          % semilato retto [km]
r = p/(1+e*cos(th));    % valore assoluto distanza dal pericentro [km]

% Vettore di stato s.r. perifocale(PF)
r_PF = r*[cos(th), sin(th), 0]';
v_PF = sqrt(mu/p)*[-sin(th), e+cos(th), 0]';

% Ruoto il vettore di stato nel sistema Geocentrico Equatoriale(ECI)

%Partenza: ECI
% Rotazione di OM intorno all'asse k
R3_OM=[cos(OM) sin(OM) 0;
      -sin(OM) cos(OM) 0;
       0       0       1];
% Rotazione di i intorno all'asse i'= N(nodo ascendente)
R1_i=[1  0      0;
      0  cos(i) sin(i);
      0 -sin(i) cos(i)];
% Rotazione di om intorno all'asse k'' = h
R3_om=[ cos(om) sin(om) 0;
       -sin(om) cos(om) 0;
        0       0       1];
% Arrivo PF

T_PF2ECI=[R3_om*R1_i*R3_OM]';     % matrice di rotazione completa PF to ECI

r = T_PF2ECI * r_PF;      %r_ECI
v = T_PF2ECI * v_PF;      %v_ECI

return




