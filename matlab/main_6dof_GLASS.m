clc; clear; close all;
%% ======================= PLOTTING (IEEE Single Column Style) ====================

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

ax_fnt  = 25;
lbl_fnt = 28;
lgd_fnt = 22;
ax_wdth = 3;
ln_wdth = 3;

%% ======================= OS4 MODEL PARAMETERS =======================
p.g  = 9.81;

p.m  = 0.24;
p.Ix = 2.3e-3;  p.Iy = 2.3e-3;  p.Iz = 4.0e-3;
p.l  = 0.20;

p.b  = 3.0e-6;
p.d  = 1.0e-7;
p.Jr = 2.0e-5;

p.omega_max = 600;
p.omega_min = 0;

%% ======================= STANDOFF CIRCLE ============================
p.xc = 0; p.yc = 0;
p.rd = 20;

%% ======================= GLASS PARAMETERS ===========================
p.kG  = 0.08;
p.dir = +1;         % +1 CCW, -1 CW

% course response limit (non-ideal)
p.kchi = 3.0;
p.wmax = 0.8;

% filter chi_d
p.tau_chid = 0.6;

%% ======================= TRANSLATION OUTER LOOP =====================
% Nominal tangential speed (this is V0 in phase-lock)
p.V0    = 2.0;
p.kv    = 2.5;
p.a_max = 8;

% radial correction term
p.krad  = 0.8;

% Tilt caps
p.phi_max = deg2rad(40);
p.th_max  = deg2rad(40);

%% ======================= PHASE-LOCK (2D reference) ==================
p.use_phase_lock = true;

% Desired angular rate for reference point.
% If you want reference to move at nominal tangential speed V0:
p.omega_d = p.V0 / p.rd;     % rad/s

% phase gain (speed scheduling); tune 0.2~1.2
p.kgamma = 0.6;

% Tangential speed saturation
p.Vmin = 0.5;
p.Vmax = 4.0;

% gamma0 will be set from initial condition so reference has SAME PHASE
p.gamma0 = 0; % placeholder; overwritten below

%% ======================= ALTITUDE LOOP ==============================
p.zd   = 10;
p.kz   = 2;
p.kvz  = 4;

%% ======================= ATTITUDE LOOP ==============================
p.kp_phi = 12; p.kd_phi = 6;
p.kp_th  = 12; p.kd_th  = 6;
p.kp_psi = 20; p.kd_psi = 10;

%% ======================= INITIAL STATE ==============================
% X = [x y z vx vy vz phi phidot theta thetadot psi psidot chi chi_d_f]
x0=40; y0=0; z0=2;
vx0=0; vy0=0; vz0=0;
phi0=0; phidot0=0;
th0 =0; thetadot0=0;
psidot0=0;

% ---- Same-phase initialization for reference + heading ----
xr = x0 - p.xc; yr = y0 - p.yc;
d0 = hypot(xr,yr);
gammaLOS0 = atan2(yr,xr);
p.gamma0 = gammaLOS0;            % <<< SAME PHASE reference at t=0

e0 = d0 - p.rd;

sigma_fun = @(e) tanh(p.kG*e);   % your style
sig0 = sigma_fun(e0);            % (no sat here; fine)

s = +1;                          % CCW
lambda0_cmd = wrapToPiLocal(s * acos(-sig0)); % consistent with your init
chi_d0 = wrapToPiLocal(gammaLOS0 + lambda0_cmd);

psi0  = chi_d0;
chi0  = chi_d0;
chidf0= chi_d0;

X0 = [x0 y0 z0 vx0 vy0 vz0 phi0 phidot0 th0 thetadot0 psi0 psidot0 chi0 chidf0]';

tspan = [0 140];
opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

[t,X] = ode45(@(t,X) f_OS4_GLASS(t,X,p), tspan, X0, opts);

%% ======================= EXTRACT ==========================
x = X(:,1); y = X(:,2); z = X(:,3);
phi = X(:,7); theta = X(:,9); psi = X(:,11);
chi = X(:,13); chi_df = X(:,14);

% Recompute desired angles for plotting (consistent with what dynamics used)
[phi_d,theta_d,psi_d,chi_d_raw] = compute_desired_history(t,X,p);

%% ======================= REFERENCE PATH (2D) ========================
gamma_ref = wrapToPiLocal(p.omega_d*t + p.gamma0);
x_ref = p.xc + p.rd*cos(gamma_ref);
y_ref = p.yc + p.rd*sin(gamma_ref);



% Compute control inputs for plotting
[U1_hist, U2_hist, U3_hist, U4_hist, Om_hist] = compute_control_history(t, X, p);

% 3D Trajectory
figure('Color','w'); ax1 = gca; hold on; grid on;
h1 = plot3(x, y, z, 'b-', 'LineWidth', ln_wdth);
plot3(x(1), y(1), z(1), 'go', 'MarkerSize', 12, 'LineWidth', ln_wdth, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
plot3(x(end), y(end), z(end), 'rs', 'MarkerSize', 12, 'LineWidth', ln_wdth, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
% Desired orbit
ang = linspace(0, 2*pi, 600);
h2 = plot3(p.xc + p.rd*cos(ang), p.yc + p.rd*sin(ang), p.zd*ones(size(ang)), 'r--', 'LineWidth', ln_wdth);
xlabel(ax1, '$x$, m', 'FontSize', lbl_fnt);
ylabel(ax1, '$y$, m', 'FontSize', lbl_fnt);
zlabel(ax1, '$z$, m', 'FontSize', lbl_fnt);
legend([h1, h2], {'UAV trajectory', 'Desired orbit'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black'; ax1.ZColor = 'black';
axis equal; box on;
view(45, 30) 

% 2D Trajectory with desired circle
figure('Color','w'); ax1 = gca; hold on; grid on;
h1 = plot(x, y, 'b-', 'LineWidth', ln_wdth);
h2 = plot(p.xc + p.rd*cos(ang), p.yc + p.rd*sin(ang), 'r--', 'LineWidth', ln_wdth);
plot(x(1), y(1), 'go', 'MarkerSize', 12, 'LineWidth', ln_wdth, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
plot(p.xc, p.yc, 'kx', 'MarkerSize', 14, 'LineWidth', ln_wdth, 'HandleVisibility', 'off');
xlabel(ax1, '$x$, m', 'FontSize', lbl_fnt);
ylabel(ax1, '$y$, m', 'FontSize', lbl_fnt);
legend([h1, h2], {'UAV trajectory', 'Desired orbit'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
axis equal; box on;

% Attitude: Roll
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, rad2deg(phi), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(phi_d), 'r--', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$\phi$, deg', 'FontSize', lbl_fnt);
legend({'$\phi$', '$\phi_d$'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

% Attitude: Pitch
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, rad2deg(theta), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(theta_d), 'r--', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$\theta$, deg', 'FontSize', lbl_fnt);
legend({'$\theta$', '$\theta_d$'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

% Attitude: Yaw
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, rad2deg(wrapToPiLocal(psi)), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(wrapToPiLocal(psi_d)), 'r--', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$\psi$, deg', 'FontSize', lbl_fnt);
legend({'$\psi$', '$\psi_d$'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

% Course commands
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, rad2deg(wrapToPiLocal(chi_d_raw)), 'g--', 'LineWidth', ln_wdth);
plot(t, rad2deg(wrapToPiLocal(chi_df)), 'r-', 'LineWidth', ln_wdth);
plot(t, rad2deg(wrapToPiLocal(chi)), 'b-', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$\chi$, deg', 'FontSize', lbl_fnt);
legend({'$\chi_d$ (raw)', '$\chi_d$ (filtered)', '$\chi$'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% ======================= CONTROL INPUT PLOTS ====================
% Thrust U1
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, U1_hist, 'b-', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$U_1$, N', 'FontSize', lbl_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

% Roll moment U2
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, U2_hist, 'r-', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$U_2$, Nm', 'FontSize', lbl_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;
ax1.XLim = [0 150];          % x from 0 to 150
ax1.YLim = [-0.1 0.1];       % y from -0.1 to 0.1

% Pitch moment U3
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, U3_hist, 'g-', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$U_3$, Nm', 'FontSize', lbl_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;
ax1.XLim = [0 150];          % x from 0 to 150
ax1.YLim = [-0.1 0.1];       % y from -0.1 to 0.1

% Yaw moment U4
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, U4_hist, 'm-', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$U_4$, Nm', 'FontSize', lbl_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;
ax1.XLim = [0 150];          % x from 0 to 150
ax1.YLim = [-0.1 0.1];       % y from -0.1 to 0.1

% All control inputs in subplots
figure('Color','w');
subplot(2,2,1); hold on; grid on;
plot(t, U1_hist, 'b-', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$U_1$, N', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;

subplot(2,2,2); hold on; grid on;
plot(t, U2_hist, 'r-', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$U_2$, Nm', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;
ax1 = gca;
ax1.XLim = [0 150];          % x from 0 to 150
ax1.YLim = [-0.1 0.1];       % y from -0.1 to 0.1

subplot(2,2,3); hold on; grid on;
plot(t, U3_hist, 'g-', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$U_3$, Nm', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;
ax1 = gca;
ax1.XLim = [0 150];          % x from 0 to 150
ax1.YLim = [-0.1 0.1]; 

subplot(2,2,4); hold on; grid on;
plot(t, U4_hist, 'm-', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$U_4$, Nm', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;
ax1 = gca;
ax1.XLim = [0 150];          % x from 0 to 150
ax1.YLim = [-0.1 0.1]; 

% Motor speeds
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, Om_hist(:,1), 'r-', 'LineWidth', ln_wdth);
plot(t, Om_hist(:,2), 'g-', 'LineWidth', ln_wdth);
plot(t, Om_hist(:,3), 'b-', 'LineWidth', ln_wdth);
plot(t, Om_hist(:,4), 'm-', 'LineWidth', ln_wdth);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$\Omega$, rad/s', 'FontSize', lbl_fnt);
legend({'$\Omega_1$', '$\Omega_2$', '$\Omega_3$', '$\Omega_4$'}, 'Location', 'best', 'FontSize', lgd_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

% Circular reference trajectory (using phase-lock reference for consistency)
gamma_ref_plot = wrapToPiLocal(p.omega_d * t + p.gamma0);
x_ref_circ = p.xc + p.rd * cos(gamma_ref_plot);
y_ref_circ = p.yc + p.rd * sin(gamma_ref_plot);

% Position Tracking (x, y, z) with circular references
figure('Color','w');
% X position
subplot(3,1,1); hold on; grid on;
plot(t, x, 'b-', 'LineWidth', ln_wdth);
plot(t, x_ref_circ, 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.8);
ylabel('$x$, m', 'FontSize', lbl_fnt*0.8);
legend({'Actual', 'Reference'}, 'Location', 'best', 'FontSize', lgd_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.8, 'LineWidth', ax_wdth*0.8);
box on;
% Y position
subplot(3,1,2); hold on; grid on;
plot(t, y, 'b-', 'LineWidth', ln_wdth);
plot(t, y_ref_circ, 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.8);
ylabel('$y$, m', 'FontSize', lbl_fnt*0.8);
legend({'Actual', 'Reference'}, 'Location', 'best', 'FontSize', lgd_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.8, 'LineWidth', ax_wdth*0.8);
box on;
% Z position
subplot(3,1,3); hold on; grid on;
plot(t, z, 'b-', 'LineWidth', ln_wdth);
yline(p.zd, 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.8);
ylabel('$z$, m', 'FontSize', lbl_fnt*0.8);
legend({'Actual', 'Reference'}, 'Location', 'best', 'FontSize', lgd_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.8, 'LineWidth', ax_wdth*0.8);
box on;
sgtitle('Position Tracking', 'FontSize', lbl_fnt);

% Attitude Tracking (φ, θ, ψ) with references
figure('Color','w');
% Roll
subplot(3,1,1); hold on; grid on;
plot(t, rad2deg(phi), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(phi_d), 'r--', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.8);
ylabel('$\phi$, deg', 'FontSize', lbl_fnt*0.8);
legend({'$\phi$', '$\phi_d$'}, 'Location', 'best', 'FontSize', lgd_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.8, 'LineWidth', ax_wdth*0.8);
box on;
% Pitch
subplot(3,1,2); hold on; grid on;
plot(t, rad2deg(theta), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(theta_d), 'r--', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.8);
ylabel('$\theta$, deg', 'FontSize', lbl_fnt*0.8);
legend({'$\theta$', '$\theta_d$'}, 'Location', 'best', 'FontSize', lgd_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.8, 'LineWidth', ax_wdth*0.8);
box on;
% Yaw
subplot(3,1,3); hold on; grid on;
plot(t, rad2deg(wrapToPiLocal(psi)), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(wrapToPiLocal(psi_d)), 'r--', 'LineWidth', ln_wdth);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.8);
ylabel('$\psi$, deg', 'FontSize', lbl_fnt*0.8);
legend({'$\psi$', '$\psi_d$'}, 'Location', 'best', 'FontSize', lgd_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.8, 'LineWidth', ax_wdth*0.8);
box on;
sgtitle('Attitude Tracking', 'FontSize', lbl_fnt);

% Combined Tracking (3 rows x 2 columns)
figure('Color','w');
% First row: X and Roll
subplot(3,2,1); hold on; grid on;
plot(t, x, 'b-', 'LineWidth', ln_wdth);
plot(t, x_ref_circ, 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$x$, m', 'FontSize', lbl_fnt*0.7);
title('Position', 'FontSize', lbl_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;

subplot(3,2,2); hold on; grid on;
plot(t, rad2deg(phi), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(phi_d), 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$\phi$, deg', 'FontSize', lbl_fnt*0.7);
title('Attitude', 'FontSize', lbl_fnt*0.8);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;

% Second row: Y and Pitch
subplot(3,2,3); hold on; grid on;
plot(t, y, 'b-', 'LineWidth', ln_wdth);
plot(t, y_ref_circ, 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$y$, m', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;

subplot(3,2,4); hold on; grid on;
plot(t, rad2deg(theta), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(theta_d), 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$\theta$, deg', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;

% Third row: Z and Yaw
subplot(3,2,5); hold on; grid on;
plot(t, z, 'b-', 'LineWidth', ln_wdth);
yline(p.zd, 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$z$, m', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;

subplot(3,2,6); hold on; grid on;
plot(t, rad2deg(wrapToPiLocal(psi)), 'b-', 'LineWidth', ln_wdth);
plot(t, rad2deg(wrapToPiLocal(psi_d)), 'r--', 'LineWidth', ln_wdth*0.7);
xlabel('$t$, s', 'FontSize', lbl_fnt*0.7);
ylabel('$\psi$, deg', 'FontSize', lbl_fnt*0.7);
set(gca, 'FontSize', ax_fnt*0.7, 'LineWidth', ax_wdth*0.7);
box on;
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
% sgtitle('Combined Position and Attitude Tracking', 'FontSize', lbl_fnt);

% Standoff error
r = hypot(X(:,1)-p.xc, X(:,2)-p.yc);
e = r - p.rd;
figure('Color','w'); ax1 = gca; hold on; grid on;
plot(t, e, 'b-', 'LineWidth', ln_wdth);
yline(0, 'k--', 'LineWidth', ln_wdth*0.7);
xlabel(ax1, '$t$, s', 'FontSize', lbl_fnt);
ylabel(ax1, '$e$, m', 'FontSize', lbl_fnt);
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

fprintf('\n=== Simulation Complete ===\n');
fprintf('Final standoff error: %.4f m\n', e(end));
fprintf('Mean standoff error:  %.4f m\n', mean(abs(e)));

%% Create 3D Trajectory GIF Animation
create_3d_trajectory_gif(x, y, z, x_ref, y_ref, p, t);

% end

%% ====================================================================
function dX = f_OS4_GLASS(t,X,p)

% States
x=X(1); y=X(2); z=X(3);
vx=X(4); vy=X(5); vz=X(6);
phi=X(7); phidot=X(8);
theta=X(9); thetadot=X(10);
psi=X(11); psidot=X(12);
chi=X(13);
chi_df=X(14);

g=p.g;

%% Geometry
dx = x - p.xc;  dy = y - p.yc;
d  = hypot(dx,dy) + 1e-9;
gamma = atan2(dy,dx);
e = d - p.rd;

%% GLASS circle: cos(lambda) = -tanh(kG e)
sigma = -tanh(p.kG*e);
sigma = max(min(sigma,1),-1);
lam_mag = acos(sigma);
lambda = (p.dir>=0)*lam_mag + (p.dir<0)*(-lam_mag);

chi_d = wrapToPiLocal(gamma + lambda);

% filter chi_d
chi_df_dot = (1/p.tau_chid) * wrapToPiLocal(chi_d - chi_df);

% non-ideal course dynamics
chi_err = wrapToPiLocal(chi_df - chi);
chi_dot_cmd = p.kchi * chi_err;
chi_dot = sat(chi_dot_cmd, p.wmax);

%% ================== Phase-locked 2D reference (same phase) ==========
gamma_ref = wrapToPiLocal(p.omega_d * t + p.gamma0);
e_gamma   = wrapToPiLocal(gamma - gamma_ref);

% Tangential and radial unit vectors
t_hat = [cos(chi); sin(chi)];
r_hat = [cos(gamma); sin(gamma)];

% Tangential speed scheduling for phase lock
Vt = p.V0;
if isfield(p,'use_phase_lock') && p.use_phase_lock
    Vt = p.V0 - p.kgamma * p.rd * e_gamma;
end
Vt = min(max(Vt, p.Vmin), p.Vmax);

% Velocity reference (tangent + radial correction)
vref = Vt * t_hat - p.krad * e * r_hat;

% velocity-to-accel
vxy = [vx; vy];
a_cmd = p.kv*(vref - vxy);
if norm(a_cmd) > p.a_max
    a_cmd = (p.a_max/norm(a_cmd))*a_cmd;
end
ax=a_cmd(1); ay=a_cmd(2);

%% accel -> desired tilt
theta_d = ( ax*cos(psi) + ay*sin(psi) )/g;
phi_d   = ( ax*sin(psi) - ay*cos(psi) )/g;
theta_d = max(min(theta_d, p.th_max),  -p.th_max);
phi_d   = max(min(phi_d,   p.phi_max), -p.phi_max);

% yaw tracks course state
psi_d = chi;

%% altitude -> thrust
z_err  = p.zd - z;
vz_err = 0 - vz;
U1_des = p.m*( g + p.kz*z_err + p.kvz*vz_err );
U1_des = max(U1_des, 0);

%% attitude -> desired ang accels -> U2,U3,U4 (one-pass gyro approx)
phi_dd_des   = p.kp_phi*(phi_d - phi)     + p.kd_phi*(0 - phidot);
theta_dd_des = p.kp_th *(theta_d-theta)   + p.kd_th *(0 - thetadot);
psi_dd_des   = p.kp_psi*(psi_d - psi)     + p.kd_psi*(0 - psidot);

Omega_est = 0;

U2_des = (p.Ix/p.l)*( phi_dd_des ...
         - thetadot*psidot*((p.Iy-p.Iz)/p.Ix) ...
         + (p.Jr/p.Ix)*thetadot*Omega_est );

U3_des = (p.Iy/p.l)*( theta_dd_des ...
         - phidot*psidot*((p.Iz-p.Ix)/p.Iy) ...
         - (p.Jr/p.Iy)*phidot*Omega_est );

U4_des = p.Iz*( psi_dd_des ...
         - phidot*thetadot*((p.Ix-p.Iy)/p.Iz) );

% rotor allocation to compute Omega for gyro term
[~, Omega_big] = allocate_rotors(U1_des, U2_des, U3_des, U4_des, p);
Omega_est = Omega_big;

% OS4 rotational accelerations
phidd = thetadot*psidot*((p.Iy-p.Iz)/p.Ix) - (p.Jr/p.Ix)*thetadot*Omega_est + (p.l/p.Ix)*U2_des;
thetadd = phidot*psidot*((p.Iz-p.Ix)/p.Iy) + (p.Jr/p.Iy)*phidot*Omega_est + (p.l/p.Iy)*U3_des;
psidd = phidot*thetadot*((p.Ix-p.Iy)/p.Iz) + (1/p.Iz)*U4_des;

%% OS4 translational accelerations
xdd = (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi))*(U1_des/p.m);
ydd = (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi))*(U1_des/p.m);
zdd = -g + (cos(phi)*cos(theta))*(U1_des/p.m);

%% pack
dX = zeros(14,1);
dX(1)=vx; dX(2)=vy; dX(3)=vz;
dX(4)=xdd; dX(5)=ydd; dX(6)=zdd;
dX(7)=phidot; dX(8)=phidd;
dX(9)=thetadot; dX(10)=thetadd;
dX(11)=psidot; dX(12)=psidd;
dX(13)=chi_dot;
dX(14)=chi_df_dot;

end

%% ====================================================================
function [Om2, Omega_big] = allocate_rotors(U1,U2,U3,U4,p)
b=p.b; d=p.d;
A = [ b,  b,  b,  b;
      0, -b,  0,  b;
     -b,  0,  b,  0;
     -d,  d, -d,  d];
rhs = [U1; U2; U3; U4];
Om2 = A\rhs;

Om2 = max(Om2, p.omega_min^2);
Om2 = min(Om2, p.omega_max^2);

Om = sqrt(Om2);
Omega_big = Om(2) + Om(4) - Om(1) - Om(3);
end

function y = sat(u, umax)
y = min(max(u,-umax),umax);
end

function ang = wrapToPiLocal(ang)
ang = mod(ang + pi, 2*pi) - pi;
end

function [U1_hist, U2_hist, U3_hist, U4_hist, Om_hist] = compute_control_history(t, X, p)
% Compute control inputs for plotting
n = numel(t);
U1_hist = zeros(n,1);
U2_hist = zeros(n,1);
U3_hist = zeros(n,1);
U4_hist = zeros(n,1);
Om_hist = zeros(n,4);

for k = 1:n
    x=X(k,1); y=X(k,2); z=X(k,3);
    vx=X(k,4); vy=X(k,5); vz=X(k,6);
    phi=X(k,7); phidot=X(k,8);
    theta=X(k,9); thetadot=X(k,10);
    psi=X(k,11); psidot=X(k,12);
    chi=X(k,13); chi_df=X(k,14);
    
    g = p.g;
    
    % Geometry
    dx = x - p.xc;  dy = y - p.yc;
    d  = sqrt(dx^2 + dy^2) + 1e-9;
    gamma = atan2(dy,dx);
    e = d - p.rd;
    
    % GLASS circle
    sigma = -tanh(p.kG*e);
    sigma = max(min(sigma,1),-1);
    lam_mag = acos(sigma);
    lambda = (p.dir>=0)*lam_mag + (p.dir<0)*(-lam_mag);
    chi_d = wrapToPiLocal(gamma + lambda);
    
    % Phase-locked reference
    gamma_ref = wrapToPiLocal(p.omega_d * t(k) + p.gamma0);
    e_gamma   = wrapToPiLocal(gamma - gamma_ref);
    
    % Tangential and radial unit vectors
    t_hat = [cos(chi); sin(chi)];
    r_hat = [cos(gamma); sin(gamma)];
    
    % Tangential speed scheduling for phase lock
    Vt = p.V0;
    if isfield(p,'use_phase_lock') && p.use_phase_lock
        Vt = p.V0 - p.kgamma * p.rd * e_gamma;
    end
    Vt = min(max(Vt, p.Vmin), p.Vmax);
    
    % Velocity reference (tangent + radial correction)
    vref = Vt * t_hat - p.krad * e * r_hat;
    
    vxy = [vx; vy];
    a_cmd = p.kv*(vref - vxy);
    if norm(a_cmd) > p.a_max
        a_cmd = (p.a_max/norm(a_cmd))*a_cmd;
    end
    ax = a_cmd(1); ay = a_cmd(2);
    
    % Desired tilt
    theta_d = (ax*cos(psi) + ay*sin(psi))/g;
    phi_d   = (ax*sin(psi) - ay*cos(psi))/g;
    theta_d = max(min(theta_d, p.th_max), -p.th_max);
    phi_d   = max(min(phi_d, p.phi_max), -p.phi_max);
    psi_d = chi;
    
    % Altitude -> U1
    z_err  = p.zd - z;
    vz_err = 0 - vz;
    U1_des = p.m*(g + p.kz*z_err + p.kvz*vz_err);
    U1_des = max(U1_des, 0);
    
    % Attitude control -> U2, U3, U4
    phi_dd_des   = p.kp_phi*(phi_d - phi) + p.kd_phi*(0 - phidot);
    theta_dd_des = p.kp_th*(theta_d - theta) + p.kd_th*(0 - thetadot);
    psi_dd_des   = p.kp_psi*(psi_d - psi) + p.kd_psi*(0 - psidot);
    
    Omega_est = 0;
    
    U2_des = (p.Ix/p.l)*(phi_dd_des - thetadot*psidot*((p.Iy-p.Iz)/p.Ix) + (p.Jr/p.Ix)*thetadot*Omega_est);
    U3_des = (p.Iy/p.l)*(theta_dd_des - phidot*psidot*((p.Iz-p.Ix)/p.Iy) - (p.Jr/p.Iy)*phidot*Omega_est);
    U4_des = p.Iz*(psi_dd_des - phidot*thetadot*((p.Ix-p.Iy)/p.Iz));
    
    % Rotor allocation
    [Om2, ~] = allocate_rotors(U1_des, U2_des, U3_des, U4_des, p);
    
    U1_hist(k) = U1_des;
    U2_hist(k) = U2_des;
    U3_hist(k) = U3_des;
    U4_hist(k) = U4_des;
    Om_hist(k,:) = sqrt(max(Om2, 0))';
end
end

function [phi_d,theta_d,psi_d,chi_d_raw] = compute_desired_history(t,X,p)
% recompute desired angles using same logic as dynamics for plotting
n = numel(t);
phi_d=zeros(n,1); theta_d=zeros(n,1); psi_d=zeros(n,1); chi_d_raw=zeros(n,1);

for k=1:n
    x=X(k,1); y=X(k,2); vx=X(k,4); vy=X(k,5);
    psi=X(k,11); chi=X(k,13); chi_df=X(k,14);

    dx=x-p.xc; dy=y-p.yc;
    d = sqrt(dx^2+dy^2)+1e-9;
    gamma = atan2(dy,dx);
    e = d - p.rd;

    sigma=-tanh(p.kG*e); sigma=max(min(sigma,1),-1);
    lam_mag=acos(sigma);
    lambda = (p.dir>=0)*lam_mag + (p.dir<0)*(-lam_mag);
    chi_d = wrapToPiLocal(gamma+lambda);
    chi_d_raw(k)=chi_d;

    % Phase-locked reference
    gamma_ref = wrapToPiLocal(p.omega_d * t(k) + p.gamma0);
    e_gamma   = wrapToPiLocal(gamma - gamma_ref);
    
    % Tangential and radial unit vectors
    t_hat=[cos(chi); sin(chi)];
    r_hat=[cos(gamma); sin(gamma)];
    
    % Tangential speed scheduling for phase lock
    Vt = p.V0;
    if isfield(p,'use_phase_lock') && p.use_phase_lock
        Vt = p.V0 - p.kgamma * p.rd * e_gamma;
    end
    Vt = min(max(Vt, p.Vmin), p.Vmax);
    
    % Velocity reference (tangent + radial correction)
    vref = Vt * t_hat - p.krad * e * r_hat;

    vxy=[vx;vy];
    a_cmd=p.kv*(vref-vxy);
    if norm(a_cmd)>p.a_max, a_cmd=(p.a_max/norm(a_cmd))*a_cmd; end
    ax=a_cmd(1); ay=a_cmd(2);

    theta_d(k)=(ax*cos(psi)+ay*sin(psi))/p.g;
    phi_d(k)  =(ax*sin(psi)-ay*cos(psi))/p.g;

    theta_d(k)=max(min(theta_d(k),p.th_max),-p.th_max);
    phi_d(k)=max(min(phi_d(k),p.phi_max),-p.phi_max);

    psi_d(k)=chi;
end
end

function create_3d_trajectory_gif(x, y, z, x_ref, y_ref, p, t)
% Create animated GIF of 3D trajectory
fprintf('\nCreating 3D trajectory animation...\n');

% Setup figure
fig = figure('Color','w', 'Position',[100 100 1200 800]);
ax = gca;
hold on; grid on; axis equal;

% Plot desired orbit (full circle)
ang = linspace(0, 2*pi, 600);
plot3(p.xc + p.rd*cos(ang), p.yc + p.rd*sin(ang), p.zd*ones(size(ang)), ...
    'r--', 'LineWidth', 2, 'DisplayName','Desired orbit');

% Plot target center
plot3(p.xc, p.yc, p.zd, 'kx', 'MarkerSize', 12, 'LineWidth', 2, ...
    'HandleVisibility', 'off');

% Initialize trajectory and markers
traj_line = plot3(NaN, NaN, NaN, 'b-', 'LineWidth', 3, 'DisplayName','UAV trajectory');
start_marker = plot3(x(1), y(1), z(1), 'go', 'MarkerSize', 12, ...
    'LineWidth', 2, 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
current_marker = plot3(NaN, NaN, NaN, 'ro', 'MarkerSize', 10, ...
    'LineWidth', 2, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

% Reference trajectory (phase-locked)
ref_line = plot3(NaN, NaN, NaN, 'm--', 'LineWidth', 1.5, 'DisplayName','Reference trajectory');
ref_marker = plot3(NaN, NaN, NaN, 'mo', 'MarkerSize', 8, ...
    'LineWidth', 1.5, 'MarkerFaceColor', 'm', 'HandleVisibility', 'off');

% Set axis properties
xlabel('$x$, m','FontSize',lbl_fnt,'Interpreter','latex');
ylabel('$y$, m','FontSize',lbl_fnt,'Interpreter','latex');
zlabel('$z$, m','FontSize',lbl_fnt,'Interpreter','latex');
title('6DOF Quadrotor GLASS Guidance Animation','FontSize',lbl_fnt,'Interpreter','latex');
lgd = legend('Location','best');
lgd.FontSize = lgd_fnt;
ax.FontSize = ax_fnt;
ax.LineWidth = ax_wdth;
ax.XColor = 'black'; ax.YColor = 'black'; ax.ZColor = 'black';
box on;

% Set view angle
view(45, 30);

% Set axis limits
margin = p.rd * 0.2;
xlim([p.xc - p.rd - margin, p.xc + p.rd + margin]);
ylim([p.yc - p.rd - margin, p.yc + p.rd + margin]);
zlim([min(z)-margin, max(z)+margin]);

% Animation parameters
gif_filename = '6dof_trajectory_animation.gif';
n_frames = min(150, length(t));  % Limit frames for reasonable file size
frame_indices = round(linspace(1, length(t), n_frames));

% Preallocate GIF frames
gif_frames = {};

for frame_idx = 1:length(frame_indices)
    idx = frame_indices(frame_idx);
    t_current = t(idx);
    
    % Update UAV trajectory
    traj_line.XData = x(1:idx);
    traj_line.YData = y(1:idx);
    traj_line.ZData = z(1:idx);
    
    % Update current position
    current_marker.XData = x(idx);
    current_marker.YData = y(idx);
    current_marker.ZData = z(idx);
    
    % Update reference trajectory
    ref_line.XData = x_ref(1:idx);
    ref_line.YData = y_ref(1:idx);
    ref_line.ZData = p.zd * ones(1, idx);
    
    % Update reference position
    ref_marker.XData = x_ref(idx);
    ref_marker.YData = y_ref(idx);
    ref_marker.ZData = p.zd;
    
    % Update time display
    time_text = sprintf('Time = %.1f s', t_current);
    if isfield(ax, 'TimeAnnotation')
        delete(ax.TimeAnnotation);
    end
    ax.TimeAnnotation = annotation('textbox', [0.02, 0.95, 0.15, 0.05], ...
        'String', time_text, 'FontSize', 14, 'BackgroundColor', 'white', ...
        'EdgeColor', 'none', 'Interpreter', 'latex');
    
    % Capture frame
    drawnow;
    frame = getframe(fig);
    gif_frames{end+1} = frame.cdata;
    
    % Progress indicator
    if mod(frame_idx, 20) == 0
        fprintf('Progress: %d/%d frames\n', frame_idx, length(frame_indices));
    end
end

% Save as GIF
fprintf('Saving GIF animation...\n');
for i = 1:length(gif_frames)
    [imind, cm] = rgb2ind(gif_frames{i}, 256);
    if i == 1
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

fprintf('Animation saved as: %s\n', gif_filename);
close(fig);
end
