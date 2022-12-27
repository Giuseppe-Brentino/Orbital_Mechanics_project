function [C,Ceq] = nonlincon(u,i_dep,i_fb,i_arr)

Ceq = [];

%% Earth data
[kep_dep, mu] = uplanet(u(1), i_dep);

[r_dep, v_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3),...
    kep_dep(4), kep_dep(5), kep_dep(6), mu);

%% Saturn Data
[kep_fb, ~] = uplanet( u(2), i_fb);

[r_fb, v_fb] = kep2car(kep_fb(1), kep_fb(2), kep_fb(3),...
    kep_fb(4), kep_fb(5), kep_fb(6), mu);

%% First heliocentric leg parameters

time_dep = datetime( mjd20002date( u(1) ) );
time_fb = datetime( mjd20002date( u(2) ) );

TOF = seconds(diff([time_dep; time_fb]));

[~,~,~,~,VI,V_before_fb,~,~] = lambertMR(r_dep, r_fb,TOF, mu, 0, 0, 2);

%%  Asteroids parameter

[kep_arr, ~] = ephNEO( u(3), i_arr);

[r_arr, v_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3),...
    kep_arr(4), kep_arr(5), kep_arr(6), mu);

%% Second heliocentric leg parameters

time_arr = datetime( mjd20002date( u(3) ) );

TOF = seconds(diff([time_fb; time_arr]));

[~,~,~,~,V_after_fb,VF,~,~] = lambertMR(r_fb, r_arr,TOF, mu, 0, 0, 2);

%% DeltaV flyby

mu_Sat = astroConstants(16);

vinf_m = V_before_fb'-v_fb;
vinf_p = V_after_fb'-v_fb;

% Turning angle
delta = acos( vinf_m'*vinf_p / ( norm(vinf_m)*norm(vinf_p) ) );

% Radius of pericentre
fun = @(rp) 1e4* ( asin( 1 / (1+( rp*norm(vinf_m)^2) / mu_Sat) ) + ...
    asin(1/(1+(rp*norm(vinf_p)^2)/mu_Sat)) - delta );

R_Sat = astroConstants(26);

opt = optimoptions('fsolve','OptimalityTolerance',1e-8);
rp = fsolve(fun, R_Sat);

C(1) = -rp + R_Sat;
C(2) = -R_Sat*906.9 + rp;

end

