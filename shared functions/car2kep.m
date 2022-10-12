function [a, e, i, OM, om, th] = car2kep(r, v, settings)
%
% PROTOTYPE
% [a, e, i, OM, om, th] = car2kep(r, v, mu)
% convertion from carthesian coordinates to keplerian elements
% 
% ------------------------------INPUT-------------------------------------
% r = [km]  position vector [1x3] 
% v = [km/s] velocity vector [1x3] 
% mu = [km^3/s^2] constante di gravitazione planetaria
%
% ------------------------------OUTPUT------------------------------------
% a =  [km]  semi-major axis  
% e =  [1x1] eccentricity       
% i =  [1x1] inclinazione       
% OM = [rad] RAAN
% om = [1x1] pericentre anomaly
% th = [rad] true anomaly

% We are using ECEI (Earth Centered Equatorial Inertial Frame), (i=gamma, k=north pole)

mu=settings.mu;

r_mod=norm(r);                                         % modulo vettori posizione e velocità   
v_mod=norm(v);

h=cross(r,v);                                          % momento angolare specifico 
h_mod=norm(h);
 
i=acos(h(3)/h_mod);                                    % inclinazione

eps=v_mod^2/2-mu/r_mod;                                % energia meccanica specifica

a=-mu/(2*eps);                                         % semiasse maggiore

if i==0                                                % if i=0, N = x_axis for convention (OM=0) 
    N=[1,0,0];   
else
    N = cross([0,0,1]',h);                             % asse nodale 
end
    N_mod = norm(N);
    N=N/N_mod;

if N(2)>=0                                             % ascensione retta del nodo ascendente
    OM=acos(N(1));
else
    OM=2*pi-acos(N(1));
end

e_vec = 1/mu * (cross(v,h) - r/r_mod);                 % vettore eccentricità
e=norm(e_vec);                                         % eccentricità

if e==0
    e_vec = N;
end

if e_vec(3)>=0                                         % anomalia del pericentro
    om=acos(dot(N,e_vec)/(N_mod*e));
else 
    om=2*pi-acos(dot(N,e_vec)/(N_mod*e));
end


v_r=dot(v,r)/r_mod;                                   % velocità radiale come discrimine per l'anomalia vera


if v_r>=0                                             % anomalia vera
    th=acos(dot(e_vec,r)/(e*r_mod));
else 
    th=2*pi-acos(dot(e_vec,r)/(e*r_mod));
end

end
 

