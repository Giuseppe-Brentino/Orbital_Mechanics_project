
function dkep = kep_pert(~, kep, settings, drag)

% [da, de, di, dOM, dom, dth]

a  = kep(1);
e  = kep(2);
i  = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);

mu = settings.mu;

p = a * (1 - e^2);
r = p / (1 + e*cos(theta));
h = sqrt(p * mu);

settings.ref_sys = 'RSW';
a_p_RSW = a_pert(settings, kep, drag);

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

