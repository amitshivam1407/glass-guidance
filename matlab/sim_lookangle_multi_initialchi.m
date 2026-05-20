function out = sim_lookangle_multi_initialchi()
%SIM_LOOKANGLE_POLAR
% Paper-consistent polar simulation of GLASS look-angle standoff guidance.
% Extended state includes explicit Cartesian kinematics (chi, x, y).
%
% ------------------------- NOTATION (PAPER) -------------------------
% Target at p_T = [xT;yT]. Relative position: rvec = p - p_T.
%   d        := ||rvec||                      (range / radial distance)
%   gammaLOS := atan2(rvec_y, rvec_x)         (LOS angle)
%   chi      := vehicle course/heading angle
%   lambda   := chi - gammaLOS                (look angle, measured from LOS)
%   e        := d - r_d                        (radial error wrt desired circle)
%
% Extended kinematics (constant speed V):
%   xdot       = V*cos(chi)
%   ydot       = V*sin(chi)
%   chidot     = omega (turn rate command)
%   e_dot      = V*cos(lambda)
%   gammaLOS_dot = (V/d)*sin(lambda)
%   lambda_dot = chi_dot - gammaLOS_dot
%
% Three heading angle cases are simulated:
%   1) Ideal: chi0 = chi_d0 (on-manifold)
%   2) +20 deg deviation from ideal
%   3) -20 deg deviation from ideal
% --------------------------------------------------------------------

clear; clc; close all;

%% ================= PAPER PLOT STYLE =================
ax_fnt  = 25;
lbl_fnt = 28;
lgd_fnt = 22;
ax_wdth = 3;
ln_wdth = 2.8;

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

%% Parameters
V   = 20;          % speed [m/s]
rd  = 200;         % desired standoff radius [m]
k   = 0.02;        % shaping gain [1/m]
Tend = 120;        % sim time [s]

% Choose sigma(e) (scalar) to match your paper's sign convention.
sigma = @(e) tanh(k*e);   % paper-consistent choice for e_dot = -V*tanh(k e)

% Orbit direction branch for lambda_c(e) = s*acos(-sigma(e)):
s = +1;  % CCW

% Inner-loop dynamics
kchi = 50.0;        % look-angle tracking gain [1/s]
wmax = 0.5;        % max turn-rate [rad/s]

%% Initial condition in Cartesian
% x0  = 450;  y0  = -250;     % UAV initial position [m]
x0  = 50;  y0  = -25;     % UAV initial position [m]
xT = 0; yT = 0;             % Target position

% Convert to polar relative to target
xr = x0 - xT;  yr = y0 - yT;
d0 = hypot(xr,yr);
gammaLOS0 = atan2(yr,xr);
e0 = d0 - rd;

% Compute ideal commanded look angle and heading
sig0 = sat(sigma(e0), 0.999999);
lambda0_cmd = wrapPi(s * acos(-sig0));  % wrap to (-pi, pi]
chi_d0 = wrapPi(gammaLOS0 + lambda0_cmd);  % ideal initial heading [radians]

%% Heading deviation cases
chi_deviations_deg = [0, +90, -60];  % degrees deviation from ideal
n_cases = length(chi_deviations_deg);
case_labels = {'$\chi_0 = \chi_{\mathrm{d0}}$ (ideal)', ...
               '$\chi_0 = \chi_{\mathrm{d0}} + 90^\circ$', ...
               '$\chi_0 = \chi_{\mathrm{d0}} - 60^\circ$'};

%% Storage for results
results = struct();
colors = lines(n_cases);

%% Settling time parameters
epsTube  = 0.5;      % epsilon standoff tube [m]

fprintf('\n===== HEADING DEVIATION ANALYSIS =====\n');
fprintf('Ideal initial heading chi_d0 = %.2f deg\n', rad2deg(chi_d0));
fprintf('Ideal initial look angle lambda0_cmd = %.2f deg\n', rad2deg(lambda0_cmd));
fprintf('%-25s %-12s %-12s %-12s %-12s\n', 'Case', 'chi0 [deg]', 'lambda0 [deg]', 'T_sim [s]', 'T_ana [s]');
fprintf('%-25s %-12s %-12s %-12s %-12s\n', '-------------------------', '------------', '-------------', '------------', '------------');

for idx = 1:n_cases
    delta_chi = deg2rad(chi_deviations_deg(idx));
    chi0 = wrapPi(chi_d0 + delta_chi);
    lambda0 = wrapPi(chi0 - gammaLOS0);
    
    % Extended state: [e; lambda; gammaLOS; chi; x; y]
    x0_state = [e0; lambda0; gammaLOS0; chi0; x0; y0];
    
    % Integrate
    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
    [t,X] = ode45(@(t,x) odefun(t,x,V,rd,sigma,s,kchi,wmax), [0 Tend], x0_state, opts);
    
    e        = X(:,1);
    lambda   = X(:,2);
    gammaLOS = X(:,3);
    chi      = X(:,4);
    x_pos    = X(:,5);
    y_pos    = X(:,6);
    d        = e + rd;
    
    % Compute derivatives for analysis
    [edot, lamdot, gamdot, chidot, xdot, ydot, omega] = deal(zeros(size(t)));
    for i=1:numel(t)
        [dx, aux] = odefun(t(i), X(i,:).', V, rd, sigma, s, kchi, wmax);
        edot(i)   = dx(1);
        lamdot(i) = dx(2);
        gamdot(i) = dx(3);
        chidot(i) = dx(4);
        xdot(i)   = dx(5);
        ydot(i)   = dx(6);
        omega(i)  = aux.omega;
    end
    
    % Analytic settling time for ideal scalar ODE
    e0_abs = abs(e0);
    if e0_abs <= epsTube
        Tana = 0;
    else
        Tana = (1/(V*k)) * log( sinh(k*e0_abs) / sinh(k*epsTube) );
    end
    
    % Simulated tube entry time
    inside = (abs(e) <= epsTube);
    idx_entry = find(inside, 1, 'first');
    if ~isempty(idx_entry)
        Tsim_entry = t(idx_entry);
    else
        Tsim_entry = NaN;
    end
    
    % Store results
    results(idx).delta_chi_deg = chi_deviations_deg(idx);
    results(idx).chi0 = chi0;
    results(idx).lambda0 = lambda0;
    results(idx).t = t;
    results(idx).x = x_pos;
    results(idx).y = y_pos;
    results(idx).e = e;
    results(idx).lambda = lambda;
    results(idx).chi = chi;
    results(idx).d = d;
    results(idx).edot = edot;
    results(idx).gamdot = gamdot;
    results(idx).lamdot = lamdot;
    results(idx).chidot = chidot;
    results(idx).xdot = xdot;
    results(idx).ydot = ydot;
    results(idx).omega = omega;
    results(idx).Tana = Tana;
    results(idx).Tsim = Tsim_entry;
    results(idx).color = colors(idx,:);
    results(idx).label = case_labels{idx};
    results(idx).x_end = x_pos(end);
    results(idx).y_end = y_pos(end);
    results(idx).e_end = e(end);
    
    fprintf('%-25s %-12.2f %-13.2f %-12.3f %-12.3f\n', case_labels{idx}, rad2deg(chi0), rad2deg(lambda0), Tsim_entry, Tana);
end

fprintf('\n');

%% Plot 1: Trajectories comparison
figure('Color','w'); ax1 = gca; hold on; axis equal; grid on;
for idx = 1:n_cases
    plot(results(idx).x, results(idx).y, 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
    % Add initial position marker (square)
    plot(results(idx).x(1), results(idx).y(1), 's', 'MarkerSize', 12, ...
        'MarkerFaceColor', results(idx).color, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
    % Add endpoint marker (circle)
    plot(results(idx).x_end, results(idx).y_end, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', results(idx).color, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end
plot(xT,yT,'kx','MarkerSize',12,'LineWidth',ln_wdth,'HandleVisibility','off');
th = linspace(0,2*pi,400);
plot(xT + rd*cos(th), yT + rd*sin(th), 'r--', 'LineWidth',ln_wdth, 'DisplayName','Desired orbit');
xlabel(ax1,'$x$, m','FontSize',lbl_fnt);
ylabel(ax1,'$y$, m','FontSize',lbl_fnt);
lgd = legend('Location','best','NumColumns',1);
lgd.FontSize = lgd_fnt;
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% Plot 2: Standoff error e(t)
figure('Color','w'); ax1 = gca; hold on; grid on;
for idx = 1:n_cases
    plot(results(idx).t, results(idx).e, 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
    % Add endpoint marker
    plot(results(idx).t(end), results(idx).e_end, 'o', 'MarkerSize', 10, ...
        'MarkerFaceColor', results(idx).color, 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end
yline(epsTube,'k--','LineWidth',ln_wdth*0.7,'DisplayName','$\epsilon$ tube');

% Add markers at epsilon tube intersection with datatips
for idx = 1:n_cases
    if ~isnan(results(idx).Tsim)
        h_marker = plot(results(idx).Tsim, epsTube, 'o', 'MarkerSize', 14, ...
            'MarkerFaceColor', results(idx).color, 'MarkerEdgeColor', 'k', ...
            'LineWidth', 2, 'HandleVisibility', 'off');
        % Add datatip to show reaching time
        dt = datatip(h_marker, results(idx).Tsim, epsTube);
        dt.FontSize = 14;
    end
end

xlabel(ax1,'$t$, s','FontSize',lbl_fnt);
ylabel(ax1,'$e$, m','FontSize',lbl_fnt);
lgd = legend('Location','best','NumColumns',1);
lgd.FontSize = lgd_fnt;
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% Plot 3: Heading angle chi(t)
figure('Color','w'); ax1 = gca; hold on; grid on;
for idx = 1:n_cases
    plot(results(idx).t, rad2deg(wrapToPi(results(idx).chi)), 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
end
xlabel(ax1,'$t$, s','FontSize',lbl_fnt);
ylabel(ax1,'$\chi$, deg','FontSize',lbl_fnt);
lgd = legend('Location','best');
lgd.FontSize = lgd_fnt;
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% Plot 4: Turn rate omega(t)
figure('Color','w'); ax1 = gca; hold on; grid on;
for idx = 1:n_cases
    plot(results(idx).t, rad2deg(results(idx).omega), 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
end
xlabel(ax1,'$t$, s','FontSize',lbl_fnt);
ylabel(ax1,'$\dot{\chi}$, deg/s','FontSize',lbl_fnt);
lgd = legend('Location','best');
lgd.FontSize = lgd_fnt;
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% Plot 5: Look angle lambda(t)
figure('Color','w'); ax1 = gca; hold on; grid on;
for idx = 1:n_cases
    plot(results(idx).t, rad2deg(wrapToPi(results(idx).lambda)), 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
end
xlabel(ax1,'$t$, s','FontSize',lbl_fnt);
ylabel(ax1,'$\lambda$, deg','FontSize',lbl_fnt);
lgd = legend('Location','best');
lgd.FontSize = lgd_fnt;
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% Plot 6: Curvature vs Standoff Error
figure('Color','w'); ax1 = gca; hold on; grid on;
for idx = 1:n_cases
    kappa = results(idx).omega / V;  % curvature [1/m]
    plot(results(idx).e, kappa, 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
end
xlabel(ax1,'$e$, m','FontSize',lbl_fnt);
ylabel(ax1,'$\kappa$, m$^{-1}$','FontSize',lbl_fnt);
lgd = legend('Location','best','NumColumns',1);
lgd.FontSize = lgd_fnt;
ax1.FontSize = ax_fnt;
ax1.LineWidth = ax_wdth;
ax1.XColor = 'black'; ax1.YColor = 'black';
box on;

%% Pack outputs
out.results = results;
out.chi_deviations_deg = chi_deviations_deg;

fprintf('Initial: d0=%.2f m, e0=%.2f m\n', d0, e0);
fprintf('Ideal initial heading chi_d0=%.4f rad (%.2f deg)\n', chi_d0, rad2deg(chi_d0));

% %% Create GIF Animation
% gif_name = 'glass_guidance.gif';
% 
% if exist(gif_name,'file')
%     delete(gif_name);
% end
% 
% % Create figure with original plotting style
% fig = figure('Color','w');
% set(fig,'Position',[100 100 1200 800]);
% ax = gca;
% 
% % Plot desired orbit
% th = linspace(0,2*pi,400);
% plot(ax, xT + rd*cos(th), yT + rd*sin(th), 'r--', 'LineWidth',ln_wdth, 'DisplayName','Desired orbit');
% plot(ax, xT, yT, 'kx', 'MarkerSize',12,'LineWidth',ln_wdth,'HandleVisibility','off');
% hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
% 
% % Use original plotting style
% xlabel(ax,'$x$, m','FontSize',lbl_fnt,'Interpreter','latex');
% ylabel(ax,'$y$, m','FontSize',lbl_fnt,'Interpreter','latex');
% title_handle = title(ax,'Time = 0.0 s','FontSize',lbl_fnt,'Interpreter','latex');
% ax.FontSize = ax_fnt;
% ax.LineWidth = ax_wdth;
% ax.XColor = 'black'; ax.YColor = 'black';
% box(ax,'on');
% 
% % Animated objects - animate all trajectories simultaneously
% traj_lines = gobjects(n_cases,1);
% uav_shapes = gobjects(n_cases,1);
% 
% % Create fixed-wing UAV shapes
% uav_size = 7;
% uav_shape_x = uav_size*[-1, 2, -1, -1];
% uav_shape_y = uav_size*[-1, 0, 1, -1];
% 
% for idx = 1:n_cases
%     uav_shapes(idx) = fill(ax, NaN, NaN, results(idx).color, 'EdgeColor', 'k', 'LineWidth', 2, ...
%         'HandleVisibility','off');
% end
% 
% legend('Location','best','FontSize',lgd_fnt);
% 
% skip = 20;
% delayTime = 0.05;
% 
% for k = 1:skip:length(results(1).x)
%     % Clear and redraw everything for instantaneous animation
%     cla(ax);
% 
%     % Draw trajectories first (will appear later in legend)
%     hold(ax,'on');
% 
%     % Redraw desired orbit last (will appear first in legend)
%     plot(ax, xT + rd*cos(th), yT + rd*sin(th), 'r--', 'LineWidth',ln_wdth, 'DisplayName','Desired orbit');
% 
%     % Draw all trajectories up to current point
%     for idx = 1:n_cases
%         % Get current trajectory data
%         x_anim = results(idx).x;
%         y_anim = results(idx).y;
% 
%         % Limit to available points
%         current_k = min(k, length(x_anim));
% 
%         % Draw trajectory up to current point
%         plot(ax, x_anim(1:current_k), y_anim(1:current_k), 'LineWidth', ln_wdth, 'Color', results(idx).color, ...
%             'DisplayName', results(idx).label);
% 
%         % Calculate heading angle for UAV orientation
%         if current_k > 1
%             heading = atan2(y_anim(current_k) - y_anim(current_k-1), x_anim(current_k) - x_anim(current_k-1));
%         else
%             heading = 0;
%         end
% 
%         % Create rotated UAV shape
%         cos_h = cos(heading);
%         sin_h = sin(heading);
%         rotated_x = cos_h * uav_shape_x - sin_h * uav_shape_y + x_anim(current_k);
%         rotated_y = sin_h * uav_shape_x + cos_h * uav_shape_y + y_anim(current_k);
% 
%         % Draw UAV shape
%         % fill(ax, rotated_x, rotated_y, results(idx).color, 'EdgeColor', 'k', 'LineWidth', 2);
%         fill(ax, rotated_x, rotated_y, results(idx).color, ...
%     'EdgeColor', 'k', 'LineWidth', 2, ...
%     'HandleVisibility','off');
%     end
% 
%     % Get current time for title
%     current_time = results(1).t(current_k);
% 
%     % Restore axis properties
%     grid(ax,'on'); axis(ax,'equal');
%     xlabel(ax,'$x$, m','FontSize',lbl_fnt,'Interpreter','latex');
%     ylabel(ax,'$y$, m','FontSize',lbl_fnt,'Interpreter','latex');
%     title(ax,sprintf('Time = %.1f s', current_time),'FontSize',lbl_fnt,'Interpreter','latex');
%     ax.FontSize = ax_fnt;
%     ax.LineWidth = ax_wdth;
%     ax.XColor = 'black'; ax.YColor = 'black';
%     box(ax,'on');
% 
%     drawnow;
% 
%     frame = getframe(fig);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,32,'nodither');
% 
%     if k == 1
%         imwrite(A,map,gif_name,'gif','LoopCount',Inf,'DelayTime',delayTime);
%     else
%         imwrite(A,map,gif_name,'gif','WriteMode','append','DelayTime',delayTime);
%     end
% end
% 
% close(fig);
% fprintf('GIF animation saved to: %s\n', gif_name);

%% Create GIF Animation
gif_name = 'glass_guidance_inside.gif';

if exist(gif_name,'file')
    delete(gif_name);
end

fig = figure('Color','w');
set(fig,'Position',[100 100 1200 800]);
ax = gca;

th = linspace(0,2*pi,400);

skip = 20;
delayTime = 0.05;

uav_size = 7;
uav_shape_x = uav_size*[-1, 2, -1, -1];
uav_shape_y = uav_size*[-1, 0, 1, -1];

% Setup figure ONCE before the animation loop
hold(ax,'on');

% Desired orbit (static - draw once)
plot(ax, xT + rd*cos(th), yT + rd*sin(th), ...
    'r--', 'LineWidth', ln_wdth, ...
    'DisplayName','Desired orbit');

% Target / center (static - draw once)
plot(ax, xT, yT, 'kx', ...
    'MarkerSize',12, ...
    'LineWidth',ln_wdth, ...
    'HandleVisibility','off');

% Create trajectory lines ONCE (will update their data in loop)
traj_lines = gobjects(n_cases,1);
for idx = 1:n_cases
    traj_lines(idx) = plot(ax, NaN, NaN, ...
        'LineWidth', ln_wdth, ...
        'Color', results(idx).color, ...
        'DisplayName', results(idx).label);
end

% Create UAV markers ONCE (will update their positions in loop)
uav_markers = gobjects(n_cases,1);
for idx = 1:n_cases
    uav_markers(idx) = patch(ax, NaN, NaN, results(idx).color, ...
        'EdgeColor','k', 'LineWidth',2, 'HandleVisibility','off');
end

% Setup axis properties ONCE
grid(ax,'on');
axis(ax,'equal');
box(ax,'on');
xlabel(ax,'$x$, m','FontSize',lbl_fnt,'Interpreter','latex');
ylabel(ax,'$y$, m','FontSize',lbl_fnt,'Interpreter','latex');
h_title = title(ax,'Time = 0.0 s','FontSize',lbl_fnt,'Interpreter','latex');
ax.FontSize = ax_fnt;
ax.LineWidth = ax_wdth;
ax.XColor = 'black';
ax.YColor = 'black';
legend(ax,'Location','best','FontSize',lgd_fnt,'Interpreter','latex','Location','northeast');

for k = 1:skip:length(results(1).x)

    % Update trajectories and UAV markers (no clearing!)
    for idx = 1:n_cases

        x_anim = results(idx).x;
        y_anim = results(idx).y;
        current_k = min(k, length(x_anim));

        % Update trajectory line data
        set(traj_lines(idx), 'XData', x_anim(1:current_k), 'YData', y_anim(1:current_k));

        % Heading for UAV marker
        if current_k > 1
            heading = atan2(y_anim(current_k) - y_anim(current_k-1), ...
                            x_anim(current_k) - x_anim(current_k-1));
        else
            heading = 0;
        end

        cos_h = cos(heading);
        sin_h = sin(heading);

        rotated_x = cos_h*uav_shape_x - sin_h*uav_shape_y + x_anim(current_k);
        rotated_y = sin_h*uav_shape_x + cos_h*uav_shape_y + y_anim(current_k);

        % Update UAV marker position
        set(uav_markers(idx), 'XData', rotated_x, 'YData', rotated_y);
    end

    % Update title with current time
    current_time = results(1).t(current_k);
    set(h_title, 'String', sprintf('Time = %.1f s', current_time));

    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,32,'nodither');

    if k == 1
        imwrite(A,map,gif_name,'gif', ...
            'LoopCount',Inf, ...
            'DelayTime',delayTime);
    else
        imwrite(A,map,gif_name,'gif', ...
            'WriteMode','append', ...
            'DelayTime',delayTime);
    end
end

close(fig);
fprintf('GIF animation saved to: %s\n', gif_name);

%% ========== ODE FUNCTION ==========
function [dx, aux] = odefun(~, x, V, rd, sigma, s, kchi, wmax)
% Extended state: x = [e; lambda; gammaLOS; chi; xpos; ypos]
e        = x(1);
lambda   = wrapPi(x(2));
gammaLOS = wrapPi(x(3));
chi      = wrapPi(x(4));
% xpos   = x(5);  % not needed for dynamics
% ypos   = x(6);  % not needed for dynamics

d = e + rd;
d = max(d, 1e-3);

% --- GLASS mapping (circle case: a=0) ---
sig = sat(sigma(e), 0.999999);
lambda_c = wrapPi(s * acos(-sig));

% --- Inner loop: track lambda_c using turn-rate chi_dot = omega ---
elam = wrapPi(lambda_c - lambda);
omega_cmd = kchi * elam;
omega = sat(omega_cmd, wmax);

% Polar kinematics
edot = V*cos(lambda);
gamdot = (V/d)*sin(lambda);

% lambda_dot = chi_dot - gammaLOS_dot
lamdot = omega - gamdot;

% Explicit Cartesian kinematics
chidot = omega;
xdot = V*cos(chi);
ydot = V*sin(chi);

dx = [edot; lamdot; gamdot; chidot; xdot; ydot];
aux.omega = omega;
aux.lambda_c = lambda_c;

%% ========== HELPER FUNCTIONS ==========

function y = sat(u, umax)
% Saturate to [-umax,umax]
y = min(max(u, -umax), umax);

function a = wrapPi(a)
% Wrap angle to (-pi, pi]
a = mod(a + pi, 2*pi) - pi;

