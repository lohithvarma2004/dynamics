clear all; clc; close all;

%% Dynamic model of HPURV
% 06/08/2018 - Updated 04/09/2025

%% Model parameters
m = 13;       % mass in kg
g = 9.81;
rho = 1021;   % density of water

a = 0.17; b = 0.24; L = 0.740;      % dimensions of robot
xG = 0; yG = 0; zG = 0;             % position of centre of gravity
xB = 0; yB = 0; zB = 0;             % position of centre of buoyancy

IxG = 0.05; IyG = 0.45; IzG = 0.44; % principal MoI about CoG
Ixy = 0; Iyz = 0; Ixz = 0;          % product of inertias
Ix = IxG + m*(yG^2 + zG^2);
Iy = IyG + m*(xG^2 + zG^2);
Iz = IzG + m*(xG^2 + yG^2);

x1 = -0.257; y1 = 0; z1 = 0;        % Distances of each fin from CoG
x2 = 0.177; y2 = 0.097; z2 = 0.0495;
x3 = 0.177; y3 = -0.097; z3 = 0.0495;

Cl = 0.92; Cd = 1.12;               % Coeff of lift and drag
S1 = 0.098; L_f1 = 0.2;             % Surface area and length of caudal fin
S2 = 0.044; L_f2 = 0.1;
S3 = 0.044; L_f3 = 0.1;

l_data = [0.2, 0.1, 0.1]';
b_data = [1, 2, 3]';
S_data = [0.098, 0.044, 0.044]';
X = [l_data, b_data, ones(length(l_data),1)];
coef = X \ S_data;
S1 = coef(1)*l_data(1) + coef(2)*b_data(1) + coef(3);
S2 = coef(1)*l_data(2) + coef(2)*b_data(2) + coef(3);
S3 = coef(1)*l_data(3) + coef(2)*b_data(3) + coef(3);

PF1max = 5;
freqmax = 2;

%% Mass matrix
Mrb = [m 0 0 0 m*zG -m*yG;
       0 m 0 -m*zG 0 m*xG;
       0 0 m m*yG -m*xG 0;
       0 -m*zG m*yG Ix -Ixy -Ixz;
       m*zG 0 -m*xG -Ixy Iy -Iyz;
       -m*yG m*xG 0 -Ixz -Iyz Iz];

Xudot = 1.3; Yvdot = -2.3; Zwdot = -5.5;
Kpdot = -0.06; Nrdot = -2.04; Mqdot = -0.86;
D = -1 * [Xudot Yvdot Zwdot Kpdot Nrdot Mqdot];
Mad = diag(D);
M = Mrb + 0*Mad;

% Control gains (PD only) with ramp adjustment
Kp_base = 10 * Mrb;  % Base proportional gain
Kd_base = 5 * Mrb;   % Base derivative gain

% Simulation settings
T = 120;
dt = 0.005;  % Reduced time step for better accuracy
t = 0:dt:T;

nu = zeros(6, length(t));
eta = zeros(6, length(t));
[etad_init, etad_dot_init, ~] = generate_trajectory(0, '8', dt);
% Ensure initial theta is not near ±90° and stabilize orientation
etad_init(5) = 0.1;  % Small initial theta
eta(:,1) = etad_init;
phi = eta(4,1); theta = eta(5,1); psi = eta(6,1);
J1 = [cos(psi)*cos(theta) -cos(phi)*sin(psi)+cos(psi)*sin(phi)*sin(theta) sin(psi)*sin(phi)+cos(phi)*cos(psi)*sin(theta);
      sin(psi)*cos(theta) cos(phi)*cos(psi)+sin(psi)*sin(phi)*sin(theta) -cos(phi)*sin(psi)+sin(phi)*cos(psi)*sin(theta);
      -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
J2 = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
      0 cos(phi) -sin(phi);
      0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
J = [J1 zeros(3,3); zeros(3,3) J2];
% Smoother initial velocity with damping
ramp_init = 0.1;  % Moderate initial ramp value
nu(:,1) = 0.1 * pinv(J) * etad_dot_init / ramp_init;  % Damped initial velocity

%% Free flow velocity
Ux = 0.4; Uy = 0.4; Uz = 0;
free_flow_vel = [Ux; Uy; Uz];

% Preallocate arrays
f = zeros(6, length(t));
PF1 = zeros(1, length(t));
PF2 = zeros(1, length(t));
PF3 = zeros(1, length(t));
BFA1 = zeros(1, length(t));
BFA2 = zeros(1, length(t));
BFA3 = zeros(1, length(t));
alpha1 = zeros(1, length(t));
alpha2 = zeros(1, length(t));
alpha3 = zeros(1, length(t));
Fx1 = zeros(1, length(t));
Fy1 = zeros(1, length(t));
Fx2 = zeros(1, length(t));
Fz2 = zeros(1, length(t));
Fx3 = zeros(1, length(t));
Fz3 = zeros(1, length(t));
V1 = zeros(1, length(t));
V2 = zeros(1, length(t));
V3 = zeros(1, length(t));
V1x = zeros(1, length(t));
V1y = zeros(1, length(t));
V2x = zeros(1, length(t));
V2z = zeros(1, length(t));
V3x = zeros(1, length(t));
V3z = zeros(1, length(t));
NV1 = zeros(1, length(t));
NV2 = zeros(1, length(t));
NV3 = zeros(1, length(t));
NV1x = zeros(1, length(t));
NV1y = zeros(1, length(t));
NV2x = zeros(1, length(t));
NV2z = zeros(1, length(t));
NV3x = zeros(1, length(t));
NV3z = zeros(1, length(t));
freq1 = zeros(1, length(t));
freq2 = zeros(1, length(t));
freq3 = zeros(1, length(t));
freq11 = zeros(1, length(t));
freq22 = zeros(1, length(t));
freq33 = zeros(1, length(t));
ga1 = zeros(1, length(t)+1);
ga2 = zeros(1, length(t)+1);
ga3 = zeros(1, length(t)+1);
e1 = zeros(6, length(t));
etad = zeros(6, length(t));
eta_dot = zeros(6, length(t));
V = zeros(length(t), 1);
alpha = zeros(length(t), 1);
b1 = zeros(length(t), 1);
b2 = zeros(length(t), 1);
b3 = zeros(length(t), 1);
tau = zeros(6, length(t));

b1(1) = 30 * pi / 180;
b2(1) = 30 * pi / 180;
b3(1) = 30 * pi / 180;

%% Choose trajectory type
trajectory_type = '8'; % Change to '8' or 'circle' to test others

%% Simulation loop
for i = 1:length(t)-1
    u = nu(1,i); v = nu(2,i); w = nu(3,i); p = nu(4,i); q = nu(5,i); r = nu(6,i);
    x = eta(1,i); y = eta(2,i); z = eta(3,i); phi = eta(4,i); theta = eta(5,i); psi = eta(6,i);
    
    J1 = [cos(psi)*cos(theta) -cos(phi)*sin(psi)+cos(psi)*sin(phi)*sin(theta) sin(psi)*sin(phi)+cos(phi)*cos(psi)*sin(theta);
          sin(psi)*cos(theta) cos(phi)*cos(psi)+sin(psi)*sin(phi)*sin(theta) -cos(phi)*sin(psi)+sin(phi)*cos(psi)*sin(theta);
          -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];
    J2 = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
          0 cos(phi) -sin(phi);
          0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
    % Check Jacobian condition number
    if cond(J) > 1e6
        J = eye(6);  % Fallback to identity if ill-conditioned
    end
    Fg_inertial = [0; 0; -m * g];       
    Fb_inertial = [0; 0; m * g];
    Fg_body = J1' * Fg_inertial;        
    Fb_body = J1' * Fb_inertial;        
    F_rest = Fg_body + Fb_body;
    tau_g = cross([xG; yG; zG], Fg_body); 
    tau_b = cross([xB; yB; zB], Fb_body);
    tau_rest = tau_g + tau_b;
    g_eta = [F_rest; tau_rest]; 

    J = [J1 zeros(3,3); zeros(3,3) J2];
    nv = J * nu(:,i);
    U = J1 * free_flow_vel;
    Ux = U(1); Uy = U(2); Uz = U(3);
    
    [etad(:,i), etad_dot(:,i), etad_ddot(:,i)] = generate_trajectory(t(i), trajectory_type, dt);
    
    crb = [0 0 0 m*(yG*q+zG*r) -m*(xG*q-w) -m*(xG*r+v);
           0 0 0 -m*(yG*p+w) m*(zG*r+xG*p) -m*(yG*r-u);
           0 0 0 -m*(zG*p-v) -m*(xG*q+u) m*(zG*p+yG*q);
           -m*(yG*q+zG*r) m*(yG*p+w) m*(zG*p-v) 0 -q*Iyz-p*Ixz+r*Iz r*Iyz+p*Ixy-q*Iy;
           m*(xG*q-w) -m*(zG*r+xG*p) m*(xG*q+u) q*Iyz+p*Ixz-r*Iz 0 -r*Ixz-q*Ixy+p*Ix;
           m*(xG*r+v) m*(yG*r-u) -m*(zG*p+yG*q) -r*Iyz-p*Ixy+q*Iy r*Ixz+q*Ixy-p*Ix 0];
    cad = [0 0 0 0 -Zwdot*w Yvdot*v;
           0 0 0 Zwdot*w 0 -Xudot*u;
           0 0 0 -Yvdot*v Xudot*u 0;
           0 -Zwdot*w Yvdot*v 0 -Nrdot*r Mqdot*q;
           Zwdot*w 0 -Xudot*u Nrdot*r 0 -Kpdot*p;
           -Yvdot*v Xudot*u 0 -Mqdot*q Kpdot*p 0];
    cn = crb + 0*cad;
    cm = cn * nu(:,i);

    dn = [0.97*abs(u) 0 0 0 0 0;
          0 5.08*abs(v) 0 0 0 0;
          0 0 3.38*abs(w) 0 0 0;
          0 0 0 0.004*abs(p) 0 0;
          0 0 0 0 0.004*abs(q) 0;
          0 0 0 0 0 0.05*abs(r)];
    dm = dn * nu(:,i);

    e1(:,i) = etad(:,i) - eta(:,i);
    e1dot = etad_dot(:,i) - J * nu(:,i);
    U_body = J1' * free_flow_vel;  % Transform free-flow velocity to body frame
    U_body_full = [U_body; zeros(3,1)];  % Extend to 6x1 with zero angular velocities
    tau_ff = M * U_body_full;           % Feedforward to compensate disturbance
    
    % Ramp-adjusted gains
    ramp = max(0.1, 1 - exp(-t(i) / 0.5));  % Faster ramp
    Kp = ramp * Kp_base;
    Kd = ramp * Kd_base;
    tau(:,i) = J' * (M * etad_ddot(:,i) + Kd * e1dot + Kp * e1(:,i)) + tau_ff;  % PD with ramped gains

    if any(isnan(tau(:,i))) || any(isinf(tau(:,i))) || any(isnan(cm)) || any(isinf(cm)) || ...
       any(isnan(dm)) || any(isinf(dm)) || any(isnan(g_eta)) || any(isinf(g_eta))
        fprintf('Step %d: Invalid input detected - tau = %e, cm = %e, dm = %e, g_eta = %e\n', i, max(abs(tau(:,i))), max(abs(cm)), max(abs(dm)), max(abs(g_eta)));
        break;
    end

    B = [sin(b1(i)) cos(b1(i)) sin(b2(i)) cos(b2(i)) sin(b3(i)) cos(b3(i));
         cos(b1(i)) -sin(b1(i)) 0 0 0 0;
         0 0 cos(b2(i)) -sin(b2(i)) cos(b3(i)) -sin(b3(i));
         -z1*cos(b1(i)) z1*sin(b1(i)) y2*cos(b2(i)) -y2*sin(b2(i)) y3*cos(b3(i)) -y3*sin(b3(i));
         z1*sin(b1(i)) z1*cos(b1(i)) (z2*sin(b2(i))-x2*cos(b2(i))) (z2*cos(b2(i))+x2*sin(b2(i))) (z3*sin(b3(i))-x3*cos(b3(i))) (z3*cos(b3(i))+x3*sin(b3(i)));
         (-x1*cos(b1(i))-y1*sin(b1(i))) -(-x1*sin(b1(i))+y1*cos(b1(i))) -y2*sin(b2(i)) -y2*cos(b2(i)) -y3*sin(b3(i)) -y3*cos(b3(i))];
    f(:,i) = B \ tau(:,i);
    f(:,i) = max(min(f(:,i), 100), -100); % Force limit ±100 N

    PF1(i) = sqrt(f(1,i)^2 + f(2,i)^2);
    PF2(i) = sqrt(f(3,i)^2 + f(4,i)^2);
    PF3(i) = sqrt(f(5,i)^2 + f(6,i)^2);

    BFA1(i) = atan2(f(2,i), 0.0001 + f(1,i));
    BFA2(i) = atan2(f(4,i), 0.0001 + f(3,i));
    BFA3(i) = atan2(f(6,i), 0.0001 + f(5,i));

    alpha1(i) = atan((f(1,i) * Cl) / (0.0001 + f(2,i) * Cd));
    alpha2(i) = atan((f(3,i) * Cl) / (0.0001 + f(4,i) * Cd));
    alpha3(i) = atan((f(5,i) * Cl) / (0.0001 + f(6,i) * Cd));

    Fx1(i) = PF1(i) * sin(BFA1(i) + b1(i));
    Fy1(i) = PF1(i) * cos(BFA1(i) + b1(i));
    Fx2(i) = PF2(i) * sin(BFA2(i) + b2(i));
    Fz2(i) = PF2(i) * cos(BFA2(i) + b2(i));
    Fx3(i) = PF3(i) * sin(BFA3(i) + b3(i));
    Fz3(i) = PF3(i) * cos(BFA3(i) + b3(i));

    V1(i) = sqrt(abs(f(2,i) / (0.0001 + 0.5 * rho * Cl * S1 * sin(2 * alpha1(i)))));
    V2(i) = sqrt(abs(f(4,i) / (0.0001 + 0.5 * rho * Cl * S2 * sin(2 * alpha2(i)))));
    V3(i) = sqrt(abs(f(6,i) / (0.0001 + 0.5 * rho * Cl * S3 * sin(2 * alpha3(i)))));

    V1x(i) = -V1(i) * sin(b1(i));
    V1y(i) = V1(i) * cos(b1(i));
    V2x(i) = -V2(i) * sin(b2(i));
    V2z(i) = V2(i) * cos(b2(i));
    V3x(i) = -V3(i) * sin(b3(i));
    V3z(i) = V3(i) * cos(b3(i));

    NV1x(i) = V1x(i) - (nu(1,i) - Ux);
    NV1y(i) = V1y(i) + (nu(2,i) - Uy);
    NV2x(i) = V2x(i) - (nu(1,i) - Ux);
    NV2z(i) = V2z(i) + (nu(3,i) - Uz);
    NV3x(i) = V3x(i) - (nu(1,i) - Ux);
    NV3z(i) = V3z(i) + (nu(3,i) - Uz);

    NV1(i) = sqrt(NV1x(i)^2 + NV1y(i)^2);
    NV2(i) = sqrt(NV2x(i)^2 + NV2z(i)^2);
    NV3(i) = sqrt(NV3x(i)^2 + NV3z(i)^2);

    freq1(i) = NV1(i) / (2 * pi * L_f1);
    freq2(i) = NV2(i) / (2 * pi * L_f2);
    freq3(i) = NV3(i) / (2 * pi * L_f3);

    freq1(i) = max(min(freq1(i), freqmax), -freqmax);
    freq2(i) = max(min(freq2(i), freqmax), -freqmax);
    freq3(i) = max(min(freq3(i), freqmax), -freqmax);

    freq11(i) = (1 - exp(-t(i))) * freq1(i);
    freq22(i) = (1 - exp(-t(i))) * freq2(i);
    freq33(i) = (1 - exp(-t(i))) * freq3(i);

    ga1(i+1) = pi/6 * sin(2 * pi * freq11(i) * t(i));
    ga2(i+1) = pi/6 * sin(2 * pi * freq22(i) * t(i));
    ga3(i+1) = pi/6 * sin(2 * pi * freq33(i) * t(i));

    % Damped fin angles
    damping = 0.9;
    b1(i+1) = damping * b1(i) + (1 - damping) * max(min(ga1(i+1) + alpha1(i), pi/4), -pi/4);
    b2(i+1) = damping * b2(i) + (1 - damping) * max(min(ga2(i+1) + alpha2(i), pi/4), -pi/4);
    b3(i+1) = damping * b3(i) + (1 - damping) * max(min(ga3(i+1) + alpha3(i), pi/4), -pi/4);

    BFA1(i) = atan2(Fx1(i), 0.0001 + Fy1(i)) - b1(i+1);
    BFA2(i) = atan2(Fx2(i), 0.0001 + Fz2(i)) - b2(i+1);
    BFA3(i) = atan2(Fx3(i), 0.0001 + Fz3(i)) - b3(i+1);

    NV1(i) = freq11(i) * (2 * pi * L_f1);
    PF1(i) = sqrt((0.5*rho*Cl*S1*sin(2*alpha1(i))*abs(NV1(i))*NV1(i))^2 + (0.5*rho*Cd*S1*cos(2*alpha1(i))*abs(NV1(i))*NV1(i))^2);
    NV2(i) = freq22(i) * (2 * pi * L_f2);
    PF2(i) = sqrt((0.5*rho*Cl*S2*sin(2*alpha2(i))*abs(NV2(i))*NV2(i))^2 + (0.5*rho*Cd*S2*cos(2*alpha2(i))*abs(NV2(i))*NV2(i))^2);
    NV3(i) = freq33(i) * (2 * pi * L_f3);
    PF3(i) = sqrt((0.5*rho*Cl*S3*sin(2*alpha3(i))*abs(NV3(i))*NV3(i))^2 + (0.5*rho*Cd*S3*cos(2*alpha3(i))*abs(NV3(i))*NV3(i))^2);

    f(1,i) = PF1(i) * cos(BFA1(i));
    f(2,i) = PF1(i) * sin(BFA1(i));
    f(3,i) = PF2(i) * cos(BFA2(i));
    f(4,i) = PF2(i) * sin(BFA2(i));
    f(5,i) = PF3(i) * cos(BFA3(i));
    f(6,i) = PF3(i) * sin(BFA3(i));

    B = [sin(b1(i+1)) cos(b1(i+1)) sin(b2(i+1)) cos(b2(i+1)) sin(b3(i+1)) cos(b3(i+1));
         cos(b1(i+1)) -sin(b1(i+1)) 0 0 0 0;
         0 0 cos(b2(i+1)) -sin(b2(i+1)) cos(b3(i+1)) -sin(b3(i+1));
         -z1*cos(b1(i+1)) z1*sin(b1(i+1)) y2*cos(b2(i+1)) -y2*sin(b2(i+1)) y3*cos(b3(i+1)) -y3*sin(b3(i+1));
         z1*sin(b1(i+1)) z1*cos(b1(i+1)) (z2*sin(b2(i+1))-x2*cos(b2(i+1))) (z2*cos(b2(i+1))+x2*sin(b2(i+1))) (z3*sin(b3(i+1))-x3*cos(b3(i+1))) (z3*cos(b3(i+1))+x3*sin(b3(i+1)));
         (-x1*cos(b1(i+1))-y1*sin(b1(i+1))) -(-x1*sin(b1(i+1))+y1*cos(b1(i+1))) -y2*sin(b2(i+1)) -y2*cos(b2(i+1)) -y3*sin(b3(i+1)) -y3*cos(b3(i+1))];
    tau(:,i) = B * f(:,i);

    % RK4 Integration with initial limits
    k1 = M \ (tau(:,i) - cm - dm - g_eta);
    k2 = M \ (tau(:,i) + 0.5 * dt * k1 - cm - dm - g_eta);
    k3 = M \ (tau(:,i) + 0.5 * dt * k2 - cm - dm - g_eta);
    k4 = M \ (tau(:,i) + dt * k3 - cm - dm - g_eta);
    nu_dot = (k1 + 2*k2 + 2*k3 + k4) / 6;
    eta_dot = J * nu(:,i);
    eta(6,i+1) = atan2(etad_dot(2,i), etad_dot(1,i));
    if i < 10  % Limit initial steps
        nu(:,i+1) = nu(:,i) + min(0.1, dt) * nu_dot;  % Small initial step
        eta(:,i+1) = eta(:,i) + min(0.1, dt) * eta_dot;
    else
        nu(:,i+1) = nu(:,i) + dt * nu_dot;
        eta(:,i+1) = eta(:,i) + dt * eta_dot;
    end
    nu(:,i+1) = max(min(nu(:,i+1), 5), -5);
    eta(:,i+1) = max(min(eta(:,i+1), 10), -10);  % Broader bounds for position

    if any(isnan(eta(:,i+1))) || any(isinf(eta(:,i+1))) || any(abs(eta(:,i+1)) > 1e5) || ...
       any(isnan(nu(:,i+1))) || any(isinf(nu(:,i+1))) || any(abs(nu(:,i+1)) > 1e5)
        fprintf('Step %d: Simulation unstable - eta = %e, nu = %e\n', i, max(abs(eta(:,i+1))), max(abs(nu(:,i+1))));
        break;
    end

    V(i) = sqrt(u^2 + v^2);

    % Reduced frequency visualization
    if mod(i, 500) == 0  % Plot every 500 steps to reduce overhead
        figure(1);
        xt1 = [0, 0.2, 0.2] - [0.3, 0.3, 0.3]; yt1 = [0, 0.1, -0.1]; zt1 = [0, 0.2, 0.2]; X1 = [xt1', yt1', zt1'];
        xt2 = xt1; yt2 = [0, 0.1, 0.1]; zt2 = [0, 0.2, -0.2]; X2 = [xt2', yt2', zt2'];
        xt3 = xt1; yt3 = [0, -0.1, -0.1]; zt3 = [0, 0.2, -0.2]; X3 = [xt3', yt3', zt3'];
        xt4 = xt1; yt4 = yt1; zt4 = [0, -0.2, -0.2]; X4 = [xt4', yt4', zt4'];
        xt5 = [0.2, 0.2, 0.5, 0.5] - [0.3, 0.3, 0.3, 0.3]; yt5 = [-0.1, 0.1, 0.1, -0.1]; zt5 = [0.2, 0.2, 0.2, 0.2]; X5 = [xt5', yt5', zt5'];
        xt6 = xt5; yt6 = [-0.1, 0.1, 0.1, -0.1]; zt6 = [-0.2, -0.2, -0.2, -0.2]; X6 = [xt6', yt6', zt6'];
        xt7 = [0.2, 0.5, 0.5, 0.2] - [0.3, 0.3, 0.3, 0.3]; yt7 = [-0.1, -0.1, -0.1, -0.1]; zt7 = [0.2, 0.2, -0.2, -0.2]; X7 = [xt7', yt7', zt7'];
        xt8 = xt7; yt8 = [0.1, 0.1, 0.1, 0.1]; zt8 = [0.2, 0.2, -0.2, -0.2]; X8 = [xt8', yt8', zt8'];
        xt9 = [0.5, 0.5, 0.6, 0.6] - [0.3, 0.3, 0.3, 0.3]; yt9 = [-0.1, 0.1, 0.075, -0.075]; zt9 = [0.2, 0.2, 0.1, 0.1]; X9 = [xt9', yt9', zt9'];
        xt10 = xt9; yt10 = [-0.1, 0.1, 0.075, -0.075]; zt10 = [-0.2, -0.2, -0.1, -0.1]; X10 = [xt10', yt10', zt10'];
        xt11 = [0.5, 0.6, 0.6, 0.5] - [0.3, 0.3, 0.3, 0.3]; yt11 = [0.1, 0.075, 0.075, 0.1]; zt11 = [0.2, 0.1, -0.1, -0.2]; X11 = [xt11', yt11', zt11'];
        xt12 = xt11; yt12 = [-0.1, -0.075, -0.075, -0.1]; zt12 = [0.2, 0.1, -0.1, -0.2]; X12 = [xt12', yt12', zt12'];
        xb = [0.6, 0.6, 0.6, 0.6] - [0.3, 0.3, 0.3, 0.3]; yb = [0.075, -0.075, -0.075, 0.075]; zb = [0.1, 0.1, -0.1, -0.1]; Xb = [xb', yb', zb'];

        xc = [0, -0.2, -0.2]; yc = [0, 0, 0]; zc = [0, 0.1, -0.1]; Xc = [xc', yc', zc'];
        Rcaudal = [cos(ga1(i)) -sin(ga1(i)) 0; sin(ga1(i)) cos(ga1(i)) 0; 0 0 1];
        Xc = (Rcaudal * Xc')'; Xc(:,1) = Xc(:,1) - [0.3, 0.3, 0.3]';

        xpl = [0.05, 0.05, -0.05, -0.05]; ypl = [-0.1, -0.3, -0.3, -0.1]; zpl = [0, 0, 0, 0]; Xpl = [xpl', ypl', zpl'];
        Rpect = [cos(ga2(i)) 0 sin(ga2(i)); 0 1 0; -sin(ga2(i)) 0 cos(ga2(i))];
        Xpl = (Rpect * Xpl')'; Xpl(:,1) = Xpl(:,1) + [0.1, 0.1, 0.1, 0.1]';

        xpr = [0.05, 0.05, -0.05, -0.05]; ypr = [0.1, 0.3, 0.3, 0.1]; zpr = [0, 0, 0, 0]; Xpr = [xpr', ypr', zpr'];
        Rpect = [cos(ga3(i)) 0 sin(ga3(i)); 0 1 0; -sin(ga3(i)) 0 cos(ga3(i))];
        Xpr = (Rpect * Xpr')'; Xpr(:,1) = Xpr(:,1) + [0.1, 0.1, 0.1, 0.1]';

        R = [cos(psi)*cos(theta) cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
             sin(psi)*cos(theta) sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
            -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];

        X1 = (R * X1')'; X2 = (R * X2')'; X3 = (R * X3')'; X4 = (R * X4')'; X5 = (R * X5')'; X6 = (R * X6')';
        X7 = (R * X7')'; X8 = (R * X8')'; X9 = (R * X9')'; X10 = (R * X10')'; X11 = (R * X11')'; X12 = (R * X12')';
        Xb = (R * Xb')'; Xc = (R * Xc')'; Xpl = (R * Xpl')'; Xpr = (R * Xpr')';

        clf;
        fill3(X1(:,1)+x, X1(:,2)+y, X1(:,3)+z, 'y', X2(:,1)+x, X2(:,2)+y, X2(:,3)+z, 'y', ...
              X3(:,1)+x, X3(:,2)+y, X3(:,3)+z, 'y', X4(:,1)+x, X4(:,2)+y, X4(:,3)+z, 'y', ...
              X5(:,1)+x, X5(:,2)+y, X5(:,3)+z, 'y', X6(:,1)+x, X6(:,2)+y, X6(:,3)+z, 'y', ...
              X7(:,1)+x, X7(:,2)+y, X7(:,3)+z, 'y', X8(:,1)+x, X8(:,2)+y, X8(:,3)+z, 'y', ...
              X9(:,1)+x, X9(:,2)+y, X9(:,3)+z, 'y', X10(:,1)+x, X10(:,2)+y, X10(:,3)+z, 'y', ...
              X11(:,1)+x, X11(:,2)+y, X11(:,3)+z, 'y', X12(:,1)+x, X12(:,2)+y, X12(:,3)+z, 'y', ...
              Xb(:,1)+x, Xb(:,2)+y, Xb(:,3)+z, 'y', Xc(:,1)+x, Xc(:,2)+y, Xc(:,3)+z, 'r', ...
              Xpl(:,1)+x, Xpl(:,2)+y, Xpl(:,3)+z, 'r', Xpr(:,1)+x, Xpr(:,2)+y, Xpr(:,3)+z, 'r');
        hold on;
        plot3(etad(1,1:i), etad(2,1:i), etad(3,1:i), 'r-', 'LineWidth', 1);
        plot3(eta(1,1:i), eta(2,1:i), eta(3,1:i), 'b-', 'LineWidth', 1);
        xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        legend('Robot', 'Desired', 'Actual');
        grid on;
        title(['Step: ', num2str(i), ' - ', trajectory_type]);
        view(45, 45);
        axis equal;
        set(gca, 'XDir', 'reverse', 'ZDir', 'reverse');
        hold off;
        drawnow;
    end
end

%% Visualization (Video)
figure('Name', 'HPURV Simulation');
video = VideoWriter([trajectory_type '_trajectory_updated.avi']);
video.FrameRate = 30;  % Reduced frame rate to improve performance
open(video);

[etad_full, ~, ~] = generate_trajectory(t, trajectory_type, dt);

for i = 1:50:length(t)-1  % Reduced frame rate for video
    phi = eta(4,i); theta = eta(5,i); psi = eta(6,i);
    R = [cos(psi)*cos(theta) cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
         sin(psi)*cos(theta) sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
        -sin(theta) cos(theta)*sin(phi) cos(theta)*cos(phi)];
    x = eta(1,i); y = eta(2,i); z = eta(3,i);

    xt1 = [0, 0.2, 0.2] - [0.3, 0.3, 0.3]; yt1 = [0, 0.1, -0.1]; zt1 = [0, 0.2, 0.2]; X1 = [xt1', yt1', zt1'];
    xt2 = xt1; yt2 = [0, 0.1, 0.1]; zt2 = [0, 0.2, -0.2]; X2 = [xt2', yt2', zt2'];
    xt3 = xt1; yt3 = [0, -0.1, -0.1]; zt3 = [0, 0.2, -0.2]; X3 = [xt3', yt3', zt3'];
    xt4 = xt1; yt4 = yt1; zt4 = [0, -0.2, -0.2]; X4 = [xt4', yt4', zt4'];
    xt5 = [0.2, 0.2, 0.5, 0.5] - [0.3, 0.3, 0.3, 0.3]; yt5 = [-0.1, 0.1, 0.1, -0.1]; zt5 = [0.2, 0.2, 0.2, 0.2]; X5 = [xt5', yt5', zt5'];
    xt6 = xt5; yt6 = [-0.1, 0.1, 0.1, -0.1]; zt6 = [-0.2, -0.2, -0.2, -0.2]; X6 = [xt6', yt6', zt6'];
    xt7 = [0.2, 0.5, 0.5, 0.2] - [0.3, 0.3, 0.3, 0.3]; yt7 = [-0.1, -0.1, -0.1, -0.1]; zt7 = [0.2, 0.2, -0.2, -0.2]; X7 = [xt7', yt7', zt7'];
    xt8 = xt7; yt8 = [0.1, 0.1, 0.1, 0.1]; zt8 = [0.2, 0.2, -0.2, -0.2]; X8 = [xt8', yt8', zt8'];
    xt9 = [0.5, 0.5, 0.6, 0.6] - [0.3, 0.3, 0.3, 0.3]; yt9 = [-0.1, 0.1, 0.075, -0.075]; zt9 = [0.2, 0.2, 0.1, 0.1]; X9 = [xt9', yt9', zt9'];
    xt10 = xt9; yt10 = [-0.1, 0.1, 0.075, -0.075]; zt10 = [-0.2, -0.2, -0.1, -0.1]; X10 = [xt10', yt10', zt10'];
    xt11 = [0.5, 0.6, 0.6, 0.5] - [0.3, 0.3, 0.3, 0.3]; yt11 = [0.1, 0.075, 0.075, 0.1]; zt11 = [0.2, 0.1, -0.1, -0.2]; X11 = [xt11', yt11', zt11'];
    xt12 = xt11; yt12 = [-0.1, -0.075, -0.075, -0.1]; zt12 = [0.2, 0.1, -0.1, -0.2]; X12 = [xt12', yt12', zt12'];
    xb = [0.6, 0.6, 0.6, 0.6] - [0.3, 0.3, 0.3, 0.3]; yb = [0.075, -0.075, -0.075, 0.075]; zb = [0.1, 0.1, -0.1, -0.1]; Xb = [xb', yb', zb'];

    xc = [0, -0.2, -0.2]; yc = [0, 0, 0]; zc = [0, 0.1, -0.1]; Xc = [xc', yc', zc'];
    Rcaudal = [cos(ga1(i)) -sin(ga1(i)) 0; sin(ga1(i)) cos(ga1(i)) 0; 0 0 1];
    Xc = (Rcaudal * Xc')'; Xc(:,1) = Xc(:,1) - [0.3, 0.3, 0.3]';

    xpl = [0.05, 0.05, -0.05, -0.05]; ypl = [-0.1, -0.3, -0.3, -0.1]; zpl = [0, 0, 0, 0]; Xpl = [xpl', ypl', zpl'];
    Rpect = [cos(ga2(i)) 0 sin(ga2(i)); 0 1 0; -sin(ga2(i)) 0 cos(ga2(i))];
    Xpl = (Rpect * Xpl')'; Xpl(:,1) = Xpl(:,1) + [0.1, 0.1, 0.1, 0.1]';

    xpr = [0.05, 0.05, -0.05, -0.05]; ypr = [0.1, 0.3, 0.3, 0.1]; zpr = [0, 0, 0, 0]; Xpr = [xpr', ypr', zpr'];
    Rpect = [cos(ga3(i)) 0 sin(ga3(i)); 0 1 0; -sin(ga3(i)) 0 cos(ga3(i))];
    Xpr = (Rpect * Xpr')'; Xpr(:,1) = Xpr(:,1) + [0.1, 0.1, 0.1, 0.1]';

    X1 = (R * X1')'; X2 = (R * X2')'; X3 = (R * X3')'; X4 = (R * X4')'; X5 = (R * X5')'; X6 = (R * X6')';
    X7 = (R * X7')'; X8 = (R * X8')'; X9 = (R * X9')'; X10 = (R * X10')'; X11 = (R * X11')'; X12 = (R * X12')';
    Xb = (R * Xb')'; Xc = (R * Xc')'; Xpl = (R * Xpl')'; Xpr = (R * Xpr')';

    clf;
    fill3(X1(:,1)+x, X1(:,2)+y, X1(:,3)+z, 'y', X2(:,1)+x, X2(:,2)+y, X2(:,3)+z, 'y', ...
          X3(:,1)+x, X3(:,2)+y, X3(:,3)+z, 'y', X4(:,1)+x, X4(:,2)+y, X4(:,3)+z, 'y', ...
          X5(:,1)+x, X5(:,2)+y, X5(:,3)+z, 'y', X6(:,1)+x, X6(:,2)+y, X6(:,3)+z, 'y', ...
          X7(:,1)+x, X7(:,2)+y, X7(:,3)+z, 'y', X8(:,1)+x, X8(:,2)+y, X8(:,3)+z, 'y', ...
          X9(:,1)+x, X9(:,2)+y, X9(:,3)+z, 'y', X10(:,1)+x, X10(:,2)+y, X10(:,3)+z, 'y', ...
          X11(:,1)+x, X11(:,2)+y, X11(:,3)+z, 'y', X12(:,1)+x, X12(:,2)+y, X12(:,3)+z, 'y', ...
          Xb(:,1)+x, Xb(:,2)+y, Xb(:,3)+z, 'y', Xc(:,1)+x, Xc(:,2)+y, Xc(:,3)+z, 'r', ...
          Xpl(:,1)+x, Xpl(:,2)+y, Xpl(:,3)+z, 'r', Xpr(:,1)+x, Xpr(:,2)+y, Xpr(:,3)+z, 'r');
    hold on;
    plot3(etad_full(1,:), etad_full(2,:), etad_full(3,:), 'r-.', 'LineWidth', 0.5);
    plot3(eta(1,1:i), eta(2,1:i), eta(3,1:i), 'b-', 'LineWidth', 1);
    plot3([-4.5, 4.5], [-0.5, 8.5], [-0.5, 4.5], 'w.');
    set(gca, 'XDir', 'reverse', 'ZDir', 'reverse');
    view(45, 45);
    grid on;
    axis equal;
    xlabel('x, [m]'); ylabel('y, [m]'); zlabel('z, [m]');
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    hold off;
    drawnow;
    frame = getframe(gcf);
    resizedFrame = imresize(frame.cdata, [420, 560]);
    writeVideo(video, resizedFrame);
end

close(video);
close(gcf);

%% Plots
figure; plot(t, V, 'LineWidth', 2); xlabel('Time, [s]'); ylabel('Absolute resultant velocity, [m/s]'); grid on;
figure; e1 = etad(:,1:end-1) - eta(:,1:end-1); plot(t(1:end-1), e1(1,:), 'r-', t(1:end-1), e1(2,:), 'b-.', t(1:end-1), e1(3,:), 'g--', 'LineWidth', 1);
legend('x_e, [m]', 'y_e, [m]', 'z_e, [m]'); xlabel('Time, [s]'); ylabel('Position error (\eta_d - \eta)'); grid on; set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
figure; plot(t(1:end-1), eta(1:3,1:end-1)'); legend('x (m)', 'y (m)', 'z (m)'); grid on; xlabel('Time, [s]'); ylabel('\eta'); set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
figure; plot(t, nu(1:3,:)'); legend('u (m/s)', 'v (m/s)', 'w (m/s)'); grid on; xlabel('Time, [s]'); ylabel('\nu'); set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
figure; plot(t, freq11, 'r-', t, freq22, 'b-', t, freq33, 'g-'); legend('caudal fin', 'pectoral fin (port)', 'pectoral fin (starboard)'); xlabel('Time, [s]'); ylabel('Frequency, [Hz]'); grid on; set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

%% Updated trajectory function with smooth start
function [etad, etad_dot, etad_ddot] = generate_trajectory(t, trajectory_type, dt)
    if isscalar(t)
        etad = zeros(6, 1); etad_dot = zeros(6, 1); etad_ddot = zeros(6, 1); t_vec = t; idx = 1;
    else
        etad = zeros(6, length(t)); etad_dot = zeros(6, length(t)); etad_ddot = zeros(6, length(t)); t_vec = t; idx = 1:length(t);
    end

    % Smooth ramp for trajectory start (over ~0.5 seconds) with minimum value
    ramp = max(0.1, 1 - exp(-t_vec / 0.5));

    switch lower(trajectory_type)
        case '8'
            a = 5; omega = 0.1;
            for i = idx
                denom = 1 + sin(omega * t_vec(i))^2;
                etad(:,i) = [a * cos(omega * t_vec(i)) / denom; a * sin(omega * t_vec(i)) * cos(omega * t_vec(i)) / denom; 2 - 2 * cos(omega * t_vec(i)); 0; 0; 0];
                denom_deriv = 2 * omega * sin(omega * t_vec(i)) * cos(omega * t_vec(i));
                etad_dot(:,i) = [-a * omega * sin(omega * t_vec(i)) / denom - a * cos(omega * t_vec(i)) * denom_deriv / denom^2;
                                 a * omega * (cos(omega * t_vec(i))^2 - sin(omega * t_vec(i))^2) / denom - a * sin(omega * t_vec(i)) * cos(omega * t_vec(i)) * denom_deriv / denom^2;
                                 2 * omega * sin(omega * t_vec(i)); 0; 0; 0];
                denom_ddot = 2 * omega^2 * (cos(omega * t_vec(i))^2 - sin(omega * t_vec(i))^2);
                etad_ddot(:,i) = [-a * omega^2 * cos(omega * t_vec(i)) / denom - 2 * a * omega * sin(omega * t_vec(i)) * (-denom_deriv) / denom^2 - a * cos(omega * t_vec(i)) * (denom_ddot / denom^2 - 2 * denom_deriv^2 / denom^3);
                                  -a * omega^2 * 2 * sin(omega * t_vec(i)) * cos(omega * t_vec(i)) / denom - a * omega * (cos(omega * t_vec(i))^2 - sin(omega * t_vec(i))^2) * (-denom_deriv) / denom^2 - a * sin(omega * t_vec(i)) * cos(omega * t_vec(i)) * (denom_ddot / denom^2 - 2 * denom_deriv^2 / denom^3);
                                  2 * omega^2 * cos(omega * t_vec(i)); 0; 0; 0];
                etad(6,i) = atan2(etad_dot(2,i), etad_dot(1,i));
            end
        case 'm'
            amplitude = 3; frequency = 0.1; vertical_offset = 2;
            for i = idx
                etad(1,i) = t_vec(i) / 20;
                etad(2,i) = amplitude * (sin(2 * frequency * t_vec(i)) + 0.3 * sin(frequency * t_vec(i) + pi/2));
                etad(3,i) = vertical_offset + 0.5 * amplitude * sin(frequency * t_vec(i));
                etad_dot(1,i) = 0.05;
                etad_dot(2,i) = amplitude * (2 * frequency * cos(2 * frequency * t_vec(i)) + 0.3 * frequency * cos(frequency * t_vec(i) + pi/2));
                etad_dot(3,i) = 0.5 * amplitude * frequency * cos(frequency * t_vec(i));
                etad_ddot(1,i) = 0;
                etad_ddot(2,i) = amplitude * (-4 * frequency^2 * sin(2 * frequency * t_vec(i)) - 0.3 * frequency^2 * sin(frequency * t_vec(i) + pi/2));
                etad_ddot(3,i) = -0.5 * amplitude * frequency^2 * sin(frequency * t_vec(i));
                if abs(etad_dot(1,i)) > eps || abs(etad_dot(2,i)) > eps
                    etad(6,i) = atan2(etad_dot(2,i), etad_dot(1,i));
                elseif i > 1
                    etad(6,i) = etad(6,i-1);
                else
                    etad(6,i) = 0;
                end
            end
        case 'circle'
            for i = idx
                etad(:,i) = [2 * sin(0.1 * t_vec(i)); 2 - 2 * cos(0.1 * t_vec(i)); 2 - 2 * cos(0.1 * t_vec(i)); (pi/6) * sin(0.1 * t_vec(i)); (-pi/4) * sin(0.1 * t_vec(i)); (pi/3) * sin(0.1 * t_vec(i))];
                etad_dot(:,i) = [0.2 * cos(0.1 * t_vec(i)); 0.2 * sin(0.1 * t_vec(i)); 0.2 * sin(0.1 * t_vec(i)); 0.1 * (pi/6) * cos(0.1 * t_vec(i)); -0.1 * (pi/4) * cos(0.1 * t_vec(i)); 0.1 * (pi/3) * cos(0.1 * t_vec(i))];
                etad_ddot(:,i) = [-0.02 * sin(0.1 * t_vec(i)); 0.02 * cos(0.1 * t_vec(i)); 0.02 * cos(0.1 * t_vec(i)); -0.01 * (pi/6) * sin(0.1 * t_vec(i)); 0.01 * (pi/4) * sin(0.1 * t_vec(i)); -0.01 * (pi/3) * sin(0.1 * t_vec(i))];
                etad(6,i) = atan2(etad_dot(2,i), etad_dot(1,i));
            end
        otherwise
            error('Invalid trajectory type. Choose ''8'', ''m'', or ''circle''.');
    end

    % Apply smooth ramp
    etad = etad .* ramp;
    etad_dot = etad_dot .* ramp;
    etad_ddot = etad_ddot .* ramp;
end