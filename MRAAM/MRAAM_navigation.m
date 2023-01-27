clear all
close all
clc

load('A_comp_cl_dp.MAT');
load('B_comp_cl_dp.MAT');
load('C_comp_cl_dp_exe.MAT');
load('D_comp_cl_dp_exe.MAT');

% Define time vector
t0 = 0; 
tf = 8.83009;

% Target Parameters
nT = 3*32.2;

% Heading error
HE_rad = -20*pi/180;
% Pursuer Gain
% Np = [4,5,6];
Np = 5.95999; 
nc = zeros(length(Np),1);

linespecs = {'b--', 'g--', 'c--','y--','k--','m--','r-'};
for ii = 1 : length(Np)
    
    Npi = Np(ii);

    % Initial Conditions
    beta_rad = 0;
    RT1      = 20000;
    RT2      = 40000;
    RM1      = 0;
    RM2      = 40000;
    VM       = 2.5*967.51969;
    VT       = 300;
    VT1      = -VT*cos(beta_rad);
    VT2      =  VT*sin(beta_rad);

    % relative positions and velocities
    RTM1 = RT1 - RM1;
    RTM2 = RT2 - RM2;
    % line of sight angle and time derivative
    lambda = atan2( RTM2, RTM1 );

    % missile lead angle
    L = asin( VT*sin( beta_rad + lambda )/VM );

    % missile velocity components
    VM1  = VM*cos(lambda + L + HE_rad);
    VM2  = VM*sin(lambda + L + HE_rad);

    % Pointers to states [beta, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]
    sel_beta = 1;
    sel_RT1  = 2;
    sel_RT2  = 3;
    sel_RM1  = 4;
    sel_RM2  = 5;
    sel_VT1  = 6;
    sel_VT2  = 7;
    sel_VM1  = 8;
    sel_VM2  = 9;

    % Initial condition vector
    X_plant = zeros(1,8);
    y0 = [beta_rad, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2, X_plant]';

    options = odeset('abstol', 1e-13, 'reltol', 1e-10);

    % Integrate nonlinear 2-D engagement situation
    [t,y] = ode45(@nlinpronav, [t0, tf], y0, options, HE_rad, Npi, A_comp_cl_dp, B_comp_cl_dp, C_comp_cl_dp_exe, D_comp_cl_dp_exe, nT);

    % Missile and target positions in inertial coordinate system
    figure(1)
    plot(y(:,sel_RM1), y(:,sel_RM2), linespecs{ii}, 'linewidth', 0.5); hold on
    xlabel('Downrange [ft]', 'fontsize', 14);
    ylabel('Crossrange [ft]', 'fontsize', 14);
    set(gca, 'fontsize', 14);
    set(gcf, 'color', 'w');

    
    % target and missile velocity magnitudes
    VT = sqrt( y(:,sel_VT1).^2 + y(:,sel_VT2).^2 );

    % relative positions and velocities
    RTM1 = y(:,sel_RT1) - y(:,sel_RM1);
    RTM2 = y(:,sel_RT2) - y(:,sel_RM2);
    VTM1 = y(:,sel_VT1) - y(:,sel_VM1);
    VTM2 = y(:,sel_VT2) - y(:,sel_VM2);

    % relative distance
    RTM = sqrt(RTM1.^2 + RTM2.^2);
    [minRTM,index] = min(RTM);

    % line of sight angle and time derivative
    lambda     = atan2( RTM2, RTM1 );
    lambda_dot = (RTM1.*VTM2 - RTM2.*VTM1)./RTM.^2;

    % closing velocity
    VC = -(RTM1.*VTM1 + RTM2.*VTM2)./RTM;

    % Compute Actuator commands
    nc = Npi*VC.*lambda_dot;
    l = size(y,1);
    Az_achieved = zeros(l,1);
  
    for jj =1:l
        Az_achieved(jj) = C_comp_cl_dp_exe(9,:)*y(jj,10:17)' + D_comp_cl_dp_exe(9,:)*nc(jj);
    end
    figure(2)
    plot(t, nc./32.2, linespecs{ii}, 'linewidth', 2); hold on 
    xlabel('Time [s]', 'fontsize', 14);
    ylabel('Acceleration [G]', 'fontsize', 14);
    title('3 G Target Maneuver','fontsize',14);
    set(gca, 'fontsize', 14, 'xlim', [0 10]);
    set(gcf, 'color', 'w');
    grid on
    
end

 figure(1)
 plot(y(:,sel_RT1), y(:,sel_RT2), 'r--', 'linewidth', 0.5);
 h = legend('N = 5.9599','Target trajectory');
 set(h,'location','best')
 grid on
 figure(2)
 o = legend('Np = 5.9599');
 set(o,'location','best');
dele_deg = (180/pi)*y(:,12);
figure(3);
plot(t,dele_deg,'linewidth',2);
title('Fin displacement','fontsize',14);
ylabel('[degrees]','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
grid on;
dele_rate_deg = (180/pi)*y(:,13);
figure(4);
plot(t,dele_rate_deg,'linewidth',2);
title('Fin displacement rate','fontsize',14);
ylabel('[degrees/sec]','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
grid on;

figure(5);
plot(t,y(:,16),'linewidth',2);
title('\alpha_{hat}','fontsize',14);
ylabel('rad','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('\alpha_{hat}');
set(h,'location','best');
grid on;
figure(6);
plot(t,y(:,10),'linewidth',2);
title('\alpha','fontsize',14);
ylabel('rad','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('\alpha');
set(h,'location','best');
grid on;

figure(7);
plot(t,y(:,17),'linewidth',2);
title('q_{hat} ','fontsize',14);
ylabel('rad/sec]','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('q_{hat}');
set(h,'location','best');
grid on;

figure(8);
plot(t,y(:,11),'linewidth',2);
title('q ','fontsize',14);
ylabel('rad/sec]','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('q');
set(h,'location','best');
grid on;
figure(9);
plot(t,y(:,15),'linewidth',2);
title('eIAz_{hat} ','fontsize',14);
ylabel('ft/sec^2','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('eIAz_{hat}');
set(h,'location','best');
grid on;
figure(10);
plot(t,y(:,14),'linewidth',2);
title('eIAz ','fontsize',14);
ylabel('ft/sec^2','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('eIAz');
set(h,'location','best');
grid on;
figure(11);
plot(y(index,sel_RM1), y(index,sel_RM2), 'o',y(index,sel_RT1), y(index,sel_RT2), '*', 'linewidth', 2); 
hold on;
a = num2str(RTM(index));
title(['Miss distance = ',a,'ft'],'fontsize',14);
h = legend('Missile Position','Target Position');
set(h,'location','best');
xlabel('Downrange [ft]', 'fontsize', 14);
ylabel('Crossrange [ft]', 'fontsize', 14);
set(gca, 'fontsize', 14);
set(gcf, 'color', 'w');
grid on

figure(12)
plot(t,Az_achieved,'linewidth',2);
title('A_{z acheived} ','fontsize',14);
ylabel('ft/sec^2','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('Az_{acheived}');
set(h,'location','best');
grid on;
figure(12)
plot(t,Az_achieved,'linewidth',2);
title('A_{z acheived} ','fontsize',14);
ylabel('ft/sec^2','fontsize', 14);
xlabel('Time [s]', 'fontsize', 14);
h = legend('Az_{acheived}');
set(h,'location','best');
grid on;
%------------------------------The End-----------------------------------%%
