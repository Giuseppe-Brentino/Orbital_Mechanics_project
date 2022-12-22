
function dkep = kep_pert(~, kep, settings, drag)

% [da, de, di, dOM, dom, dth]

a  = kep(1);
e  = kep(2);
i  = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);

mu_E = settings.mu;

p = a * (1 - e^2);
r = p / (1 + e*cos(theta));
h = sqrt(p * mu_E);
v = sqrt( 2*mu_E/r - mu_E/a);

[r_vec, v_vec] = kep2car(a, e, i, OM, om, theta, mu_E);
a_p = a_pert(settings, r_vec, v_vec, drag);

t = v_vec/norm(v_vec);
h_vec = cross(r_vec, v_vec);
h_vers = h_vec/norm(h_vec);
n = cross(h_vers, t)/norm(cross(h_vers, t));

% rotation from Cartesian to THN
A_car2THN = [t n h_vers]';
a_p_THN = A_car2THN * a_p;

% rotation from THN to RSW
A_THN2RSW = h/p/v * [e*sin(theta),     -(1 + e*cos(theta)),     0 ;
    1 + e*cos(theta),     e*sin(theta),          0 ;
    0,                      0,                  1];
a_p_RSW = A_THN2RSW * a_p_THN;

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

