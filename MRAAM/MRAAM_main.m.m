clear all
close all
clc

t_sec = linspace(0, 10, 500);
w_rps = logspace(-1, 3, 500);
w_rps1 = logspace(-1,3,30);
% Given Aero Dynamic data
h = 40000;
Mach = 2.5;
V0 = Mach * 967.51969; % Trimmed air velocity in ft/sec
% Coupled System Dynamics
A_sys   = [-1.57      0           0           1        0;
           0          -0.5        0.17        0        1;
           -21.13     -2876.70    -2.1        -0.14    -0.05;
           -82.92     -11.22      -0.01       -0.57    0;
           -0.19      -11.86      -0.01       0       -0.57];

B_sys   = [ 0          -0.1        0;
            -0.07      0           0.11;
            -1234.7    -30.49      -1803.2;
            -4.82      -119.65     -7;
            14.84      0.27        -150.58];

% State Vector X = [alp,beta, p,q,r]
% Input vector u = [del_a,del_e,del_r];

% computing eigen values of the Pursuer
lam_sys = eig(A_sys);

% Aero data associated with short period dynamics
V_fps = 2.5*967.51969;
Zalpha_V = -1.57;
Zq_V_p1  = 1.0;
Malpha   = -82.92;
Mq       = -0.57;

Zdele_V  = -0.1;
Mdele    = -119.65;

Zalpha = Zalpha_V*V_fps;
Zdele  = Zdele_V*V_fps;

% Pointers to long. dynamics output
sel_output_alpha = 1;
sel_output_q     = 2;
sel_output_Az    = 3;

% Short Period approximation
Ap_sp   = [Zalpha_V Zq_V_p1; Malpha Mq]; nAp = size(Ap_sp, 1);
Bp_sp   = [Zdele_V; Mdele];
Cp_sp   = [1 0; 
           0 1;
           Zalpha 0]; nCp = size(Cp_sp,1);
Dp_sp   = [0; 0; Zdele];

% Output measurement data
Cp_alpha  = Cp_sp(sel_output_alpha, :);
Dp_alpha  = Dp_sp(sel_output_alpha, :);
Cp_q      = Cp_sp(sel_output_q, :);
Dp_q      = Dp_sp(sel_output_q, :);
Cp_Az     = Cp_sp(sel_output_Az, :);
Dp_Az     = Dp_sp(sel_output_Az, :);
    
% Regulated output: Az
Cp_reg = Cp_Az; nCp_reg = 1;
Dp_reg = Dp_Az; 

% Define plant as state space system
sel_outputs = [sel_output_alpha sel_output_q, sel_output_Az];
sys_plant = ss(Ap_sp, Bp_sp, Cp_sp(sel_outputs,:), ...
    Dp_sp(sel_outputs,:));

set(sys_plant, 'StateName',  {'alpha [rad]', 'q [r/s]'}, ...
               'OutputName', {'alpha [rad]', 'q [r/s]', 'Az [ft/s2]'}, ...
               'InputName',  {'dele [rad]'});

if max(real(eig(Ap_sp))) > 0
    ttd_sec = log(2)/max(real(eig(Ap_sp)));
    fprintf('Unstable, time to double = %2.2f [s]\n',ttd_sec);  
end
  
%% Second order actuator model. States: delta_e, delta_e_dot. Inputs: 
% delta_e_command, Outputs: delta_e
%--------------------------------------------------------------------------
wa = 2*pi*35;   % Natural frequency, rad/s
za = 0.71;     % Damping ratio    , n/d

Ap_act = [0, 1; -wa^2, -2*za*wa]; 
Bp_act = [0;  wa^2];       
Cp_act = [1, 0];
Dp_act = 0;

sys_act = ss(Ap_act, Bp_act, Cp_act, Dp_act);
set(sys_act, 'StateName',  {'dele [rad]', 'ddele [r/s]'}, ...
             'OutputName', {'dele [rad]'}, ...
             'InputName',  {'\delta_e^{cmd} [rad]'});

% Build servomechanism design model
%--------------------------------------------------------------------------
Aw = [zeros(nCp_reg,nCp_reg), Cp_reg;
      zeros(nAp,nCp_reg)    , Ap_sp ]; nAw = size(Aw,1);
Bw = [Dp_reg; Bp_sp]; [nBw, mBw] = size(Bw);

% Check extended system controlability
CO = ctrb(Aw, Bw);
svd_CO = svd(CO);
very_small_svd_CO = find(svd_CO <= 1e-06, 1);

if isempty(very_small_svd_CO)
    disp('Extended System Controllable');
else
    disp('Extended System NOT Controllable');
end

% State feedback of err_{Az} and alpha, and q
Cw_meas = eye(3);
Dw_meas = 0*Cw_meas*Bw;
        
sys_w = ss(Aw, Bw, Cw_meas, Dw_meas);
set(sys_w, 'StateName',  {'eIAz', 'alpha [rad]', 'q [r/s]'}, ...
           'OutputName', {'eIAz', 'alpha [rad]', 'q [r/s]'}, ...
           'InputName',  {'\delta_e^ [rad]'});

%% LQR Design loop
%--------------------------------------------------------------------------

% Specify control penalty
R = 1;

% Define penalty vector (the max and min q_eI values were determined
% through trial and error)
q_eI = logspace(-8, -4, 30); 
q_eI_lim = length(q_eI);

% Preallocate
design_data = zeros(q_eI_lim, 16);
lqr_gain = zeros(q_eI_lim,3);
closed_loop_eigs_all = zeros(5, q_eI_lim);
min_closed_loop_eigs_sysw = zeros(q_eI_lim,1);

% Pointers to closed loop system 
sel_state_alpha = 1;
sel_state_q = 2;
sel_state_dele = 3;
sel_state_ddele = 4;
sel_state_eIAz = 5;

% Plant with actuator model (analysis model)
%----------------------------------------------------------------------
sys_plant_act = sys_plant*sys_act;

Ap_sp_act = sys_plant_act.a;
Bp_sp_act = sys_plant_act.b;
Cp_sp_act = sys_plant_act.c; 
Dp_sp_act = sys_plant_act.d;

% Step options
options = stepDataOptions('StepAmplitude', 32.2);

% Get the transmission zeros for asymptotics
h_ss        = ss(Aw, Bw, [1 0 0], 0);
h_tf        = tf(h_ss);
h_tzero     = [-42.6603,42.0903]; 
% State feedback design loop
for ii = 1 : q_eI_lim
    
    % LQR PI design using all states
    Qw_eI    = q_eI(ii); 
    Qw_alpha = 0;
    Qw_q     = 0;

    Qw = diag([Qw_eI, Qw_alpha, Qw_q]);
    Rw = R; 
    
    % Check extended system observability
    OM = obsv(Aw, Qw);
    svd_OM = svd(OM);
    very_small_svd_OM = find(svd_OM <= 1e-09, 1);

    if isempty(very_small_svd_OM)
    else
        disp('(Aw,Q^(1/2)) not observable');
        return
    end

    % Solve ARE, compute gains (note design model does not include actuator
    % model)
    [Kc, Pw, ~] = lqr(Aw, Bw, Qw, Rw);
    lqr_gain(ii,:) = Kc ;


    % LQR solution check
    Ric_Residual = Pw*Aw + Aw'*Pw + Qw - Pw*Bw*Rw^(-1)*Bw'*Pw;
    norm_residual = norm(Ric_Residual);

    if norm_residual > 1e-6
        disp(['Warning: LQR ARE Solution Tolerance is: ', ...
            num2str(norm_residual)]);
    end
    
    % Dynamic controller (intergral error states, control output, plant 
    % state input )
    %----------------------------------------------------------------------
    Ac  = zeros(nCp_reg);
    Bc1 = zeros(nCp_reg,3); Bc1(1,sel_output_Az) = 1; 
    Bc2 = -eye(nCp_reg);
    Cc  = -Kc(nCp_reg);
    Dc1 = [-Kc(2:3) 0]; 
    Dc2 = zeros(nCp_reg);   
    
    % Closed loop between plant and controller (see Chapter 1.5 in Lavretsky
    % and Wise)
    %----------------------------------------------------------------------
    [Acl, Bcl, Ccl, Dcl, Bcl2, Dcl2] = SSR2Y(Ap_sp_act, Bp_sp_act, ...
        Cp_sp_act, Dp_sp_act, Ac, Bc1, Bc2, Cc, Dc1, Dc2);
    
    % Time domain analysis of closed loop control system with actuators
    %----------------------------------------------------------------------
    
    % Store the closed loop eigs and damping coefficient
    [~, Z, P] = damp(Acl);
    closed_loop_eigs_all(:,ii) = sort(P);
    min_damp_Acl = min(Z);

    % Make closed loop system a state-space object
    closed_loop_system = ss(Acl, Bcl, Ccl, Dcl);
    clp_sys_tf  = tf(closed_loop_system);
 
    % Compute step response of closed-loop system (note: step.m can
    % determine OS, US, and rise time with stepinfo).
    [y, ~, x_cl] = step(closed_loop_system, t_sec);

    % Get Az response to Az step command
    y_reg = y(:, sel_output_Az); 
     
    % Output regulation error
    alpha_err = abs( ones(size(y_reg)) - y_reg );

    % Get final error value
    fv = alpha_err(end); 

    if fv < 0.05
        % Compute 63% rise time
        e_n  = alpha_err - fv*ones(size(alpha_err)) - ...
            0.36*ones(size(alpha_err)); 
        e_n1 = abs(e_n) + e_n;
        rise_time_alpha_sec = crosst(e_n1,t_sec);   
    else       
        rise_time_alpha_sec = 0;      
    end

    % Compute 95% settling time
    e_n  = alpha_err - fv*ones(size(alpha_err)) - ...
        0.05*ones(size(alpha_err));
    e_n1 = abs(e_n) + e_n;
    settling_time_alpha_sec = crosst(e_n1, t_sec);

    % Percent undershoot
    percent_undershoot = abs(min(y_reg))*100;  

    % Percent overshoot
    percent_overshoot = ( abs(max(y_reg))-1 )*100; 

    % Max elevator displacement ( actuator displacement - fin disp)
    % x_cl(:,2)
    max_elevator_angle_deg = 180/pi*max(abs(x_cl(:, sel_state_dele)));

    % Max elevator rate ( actuator rate - fin displ rate) x_cl(:,3)
    max_elevator_rate_dps = 180/pi*max(abs(x_cl(:, sel_state_ddele)));

    wsr     = logspace(-1,3,100);
    % Complimentary-Sensitivity and Sensitivity of closed loop system with 
    % the controller in loop
    max_T_Az        = max(squeeze(mag2db(abs(freqresp(clp_sys_tf(3),w_rps)))));
    max_S_Az        = max(squeeze(mag2db(abs(freqresp(1-clp_sys_tf(3),w_rps)))));

    max_T_q         =max(squeeze(mag2db(abs(freqresp(clp_sys_tf(2),w_rps)))));
    max_S_q        = max(squeeze(mag2db(abs(freqresp(1-clp_sys_tf(2),w_rps)))));


    % Plot step response data
    figure(2)
    plot(t_sec, y_reg,'linewidth',.5); hold on 
    xlabel('Time [sec]','fontsize',14);
    ylabel('Az [ft/s]','fontsize',14);
    title('Az Closed Loop Step Response','fontsize',14);
    set(gca,'fontsize',14);
    set(gcf,'color','w');
    grid on 
 
    % Frequency domain robustness analysis
    %----------------------------------------------------------------------

    % Build loop gain at the plant input
    controller = ss(Ac, Bc1, Cc, Dc1);
    Lu      = controller*sys_plant_act; [nLu, mLu] = size(Lu.a);
    Lu.c    = -Lu.c; 
    Lu.d    = -Lu.d;
    lu_tf   =   tf(Lu);
    % Evaluate loop gain at plant input over frequency range
    Lu_eval = zeros(size(w_rps)); 
    eye_nLu = eye(nLu);
    for i = 1 :length(w_rps)
        resolvent = (1i*w_rps(i)*eye_nLu - Lu.a)\eye_nLu;
        Lu_eval(i) = Lu.d + Lu.c*resolvent*Lu.b;
    end

    % Get stability margins
    [GM, PM, VM] = stability_margins( w_rps, Lu_eval, 'off');
    GM_dB  = 20*log10(GM);
    GM_dB  = min(GM_dB(GM_dB > 0));
    PM_deg = min(PM(PM > 0));
     


     cl_eigs = squeeze(min(real(eig(Aw - Bw*Kc))));
     min_closed_loop_eigs_sysw(ii) = cl_eigs;

    % Compute crossover frequency of loop gain at plant input
    Lu_dB = 20*log10(abs(Lu_eval));
    wc_rps = crosst(Lu_dB, w_rps);
    wc_Hz = wc_rps/6.28;
     
    % Compute return difference magnitude
    RD_mag = abs(1 + Lu_eval);
    RD_dB = 20*log10(RD_mag);
    RD_min = min(RD_mag);
    RD_min_dB = 20*log10(RD_min);
     
    design_data(ii,:) = [rise_time_alpha_sec, settling_time_alpha_sec,...
       percent_overshoot, percent_undershoot, max_elevator_angle_deg, ...
       max_elevator_rate_dps, wc_Hz, GM_dB, PM_deg, RD_min, RD_min_dB, ...
       min_damp_Acl,max_S_Az,max_T_Az,max_S_q,max_T_q];

    fprintf('Design Iteration %2.0f complete.\n',ii);

    pause(0.01)
    
end

set(closed_loop_system, ...
    'StateName',  {'alpha [rad]', 'q [r/s]', 'dele [rad]', 'ddele [r/s]', ...
                    'eIalpha [r-sec]'}, ...
    'OutputName', {'alpha [rad]', 'q [r/s]', 'Az [ft/s2]'}, ...
    'InputName',  {'alpha_{cmd} [rad]'})
    
%% Design charts
%--------------------------------------------------------------------------
sel_dd_taur = 1;
sel_dd_taus = 2;
sel_dd_OS   = 3;
sel_dd_US   = 4;
sel_dd_dele = 5;
sel_dd_ddele= 6;
sel_dd_wc   = 7;
sel_dd_GM   = 8;
sel_dd_PM   = 9;
sel_dd_RDmin= 10;
sel_dd_VM   = 11;
sel_dd_damp = 12;
sel_dd_S_Az = 13;
sel_dd_T_Az = 14;
sel_dd_S_q  = 15;
sel_dd_T_q  = 16;

figure(3);
plot(design_data(:,sel_dd_wc), design_data(:,sel_dd_OS),'*','linewidth',0.75); hold on
plot(design_data(:,sel_dd_wc), design_data(:,sel_dd_US),'r*','linewidth',0.75); hold on 
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel('[%]','fontsize', 14);
title('Design Chart: Percent Overshoot and Undershoot','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 
h = legend('%OS','%US');
set(h,'location','best');

figure(4);
plot(design_data(:,sel_dd_wc), design_data(:,sel_dd_dele)*644,'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel({'fin Displacement (Deg) ','[(20g ft/s^2 Az)]'},'fontsize', 14);
title('Design Chart: Fin displacement','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(5);
plot(design_data(:,sel_dd_wc), design_data(:,sel_dd_ddele)*644,'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel({'Fin Displacement Rate ( Deg/sec)','[(20g ft/s^2 Az)]'},'fontsize', 14);
title('Design Chart: Fin Displacement Rate','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(6);
plot(design_data(:,sel_dd_wc), design_data(:,sel_dd_GM),'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel('Gain Margin [dB]','fontsize', 14);
title('Design Chart: Gain Margin at Plant Input','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(7);
plot(design_data(:,sel_dd_wc), design_data(:,sel_dd_PM),'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel('Phase Margin [deg]','fontsize', 14);
title('Design Chart: Phase Margin at Plant Input','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(8)
semilogx(w_rps1, design_data(:, sel_dd_S_Az),'b*','linewidth',0.75); hold on
semilogx(w_rps1,design_data(:,sel_dd_T_Az),'r*','linewidth',0.75);
xlabel('Frequency, \omega rad/sec', 'fontsize', 14);
ylabel('|S| & |T| [dB]', 'fontsize', 14);
title('Sensitivity & Co-sensitivity of Az ','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on

figure(9);
semilogx(w_rps1, design_data(:, sel_dd_S_q),'b*','linewidth',0.75); hold on
semilogx(w_rps1, design_data(:, sel_dd_T_q),'r*','linewidth',0.75); hold on
xlabel('Frequency, \omega rad/sec', 'fontsize', 14);
ylabel('|S| & |T| of q [dB]' , 'fontsize', 14);
title('Sensitivity and co-sensitivity q ','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on

figure(11)
semilogx(q_eI(1,:), design_data(:,sel_dd_wc)*6.28,'b*','linewidth',0.75); hold on
xlabel('Penalty q_ei', 'fontsize', 14);
ylabel('CrossOver Frequency','fontsize', 14);
title('Design Chart: Crossover frequency','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

%% Based on design chart analysis, design point selected is
dp_select = 19;


figure(3); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_OS),'ro','markersize',8,'linewidth',0.75);
           plot(design_data(dp_select,sel_dd_wc), ...
           design_data(dp_select, sel_dd_US),'ro','markersize',8,'linewidth',0.75);
figure(4); plot(design_data(dp_select,sel_dd_wc), ...
644*design_data(dp_select, sel_dd_dele),'ro','markersize',8,'linewidth',0.75);
figure(5); plot(design_data(dp_select,sel_dd_wc), ...
644*design_data(dp_select, sel_dd_ddele),'ro','markersize',8,'linewidth',0.75);
figure(6); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_GM),'ro','markersize',8,'linewidth',0.75);
figure(7); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_PM),'ro','markersize',8,'linewidth',0.75);
figure(8); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_S_Az),'ro','markersize',8,'linewidth',0.75);
figure(9); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_T_Az),'ro','markersize',8,'linewidth',0.75);
figure(10); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_S_q),'ro','markersize',8,'linewidth',0.75);
figure(11); plot(design_data(dp_select,sel_dd_wc), ...
design_data(dp_select, sel_dd_T_q),'ro','markersize',8,'linewidth',0.75);

%% Select design point from charts and recreate design
%--------------------------------------------------------------------------

% close all design figures
%%close all

% LQR PI design using all states
Qw_eI    = q_eI(dp_select);  
Qw_alpha = 0;
Qw_q     = 0;

Qw = diag([Qw_eI, Qw_alpha, Qw_q]);
Rw = R;  

% Solve ARE, compute gains (note design model does not include actuator
% model)
[Kc, Pw, ~] = lqr(Aw, Bw, Qw, Rw);

% LQR solution check
Ric_Residual = Pw*Aw + Aw'*Pw + Qw - Pw*Bw*Rw^(-1)*Bw'*Pw;
norm_residual = norm(Ric_Residual);

if norm_residual > 1e-6
    disp(['Warning: LQR ARE Solution Tolerance is: ', ...
        num2str(norm_residual)]);
end

% Dynamic controller (intergral error states, control output, plant 
% state input )
%----------------------------------------------------------------------
Ac = zeros(nCp_reg);
Bc1 = zeros(nCp_reg,3); Bc1(1,sel_output_Az) = 1; 
Bc2 = -eye(nCp_reg);
Cc = -Kc(nCp_reg);
Dc1 = [-Kc(2:3) 0]; 
Dc2 = zeros(nCp_reg);  

% Closed loop between plant and controller (see Chapter 1.5 in Lavretsky
% and Wise)
%----------------------------------------------------------------------
[Acl, Bcl, Ccl, Dcl, Bcl2, Dcl2] = SSR2Y(Ap_sp_act, Bp_sp_act, ...
    Cp_sp_act, Dp_sp_act, Ac, Bc1, Bc2, Cc, Dc1, Dc2);

% Time domain analysis of closed loop control system with actuators
%----------------------------------------------------------------------


% Make closed loop system a state-space object
closed_loop_system = ss(Acl, Bcl*32.2, Ccl, Dcl*32.2); 
% closed_loop_system = ss(Acl, Bcl*644, Ccl, Dcl*644); 

sel_state_dele    = 3;
sel_state_deledot = 4;

% Compute step response of closed-loop system
t_sec = linspace(0,10,2000);
[y, ~, x_cl] = step(closed_loop_system, t_sec);

% Get AoA response to AoA step command
y_reg = y(:, sel_output_Az); 

% Output regulation error
alpha_err = abs( ones(size(y_reg)) - y_reg );

% Get final error value
fv = alpha_err(end); 

if fv < 0.05
    % Compute 63% rise time
    e_n  = alpha_err - fv*ones(size(alpha_err)) - ...
        0.36*ones(size(alpha_err)); 
    e_n1 = abs(e_n) + e_n;
    rise_time_alpha_sec = crosst(e_n1,t_sec);   
else       
    rise_time_alpha_sec = 0;      
end

% Compute 95% settling time
e_n  = alpha_err - fv*ones(size(alpha_err)) - ...
    0.05*ones(size(alpha_err));
e_n1 = abs(e_n) + e_n;
settling_time_alpha_sec = crosst(e_n1, t_sec);

% Build loop gain at the plant input
controller = ss(Ac, Bc1, Cc, Dc1);
set(controller,'StateName','eI_alpha');
set(controller,'InputName',{'alpha [rad]','q [r/s]', 'Az [ft/s2]'});
set(controller,'OutputName','\delta_e^{cmd} [rad]');
Lu = controller*sys_plant_act; [nLu, mLu] = size(Lu.a);
Lu.c = -Lu.c; 
Lu.d = -Lu.d; 

% Get stability margins
[GM, PM, VM] = stability_margins( w_rps, Lu_eval, 'off');
GM_dB  = 20*log10(GM);
GM_dB  = min(GM_dB(GM_dB > 0));
PM_deg = min(PM(PM > 0));

% Evaluate loop gain at plant input over frequency range
Lu_eval = zeros(size(w_rps)); 
eye_nLu = eye(nLu);
for i = 1 :length(w_rps)
    resolvent = (1i*w_rps(i)*eye_nLu - Lu.a)\eye_nLu;
    Lu_eval(i) = Lu.d + Lu.c*resolvent*Lu.b;
end
%% Bode diagram
figure(13)
margin(Lu)
set(gcf,'color','w');


%% Part 3 : Leuenberger Observer
% a) Determine the observer gains with ARE
% b) Analyze the frequency domain properties with the observer in loop
% c) Build closed loop system and evaluate time domain properties

% Cw_meas         =   [0 0 1; 1 0 0];
Cw_meas         =   [1 0 0; 0 0 1];

nCw_meas        =   size(Cw_meas,1);
K_lqr = lqr_gain(dp_select,:);
% Observer Penalty matrices
q_obs = logspace(0,4,30);
q_obs_lim = length(q_obs);
% wc_design = design_data(dp_select,7)*6.28;
wc_design = design_data(dp_select,7)*6.28;

% Memory allocation
design_data_obs = zeros(q_obs_lim,12);
min_closed_loop_eigs_sys_comp = zeros(q_obs_lim,1);
for jj = 1:q_obs_lim
    w  = q_obs(jj)*eye(3);
    v  = 1*[1 0;0 0.1];

% Define state space integrator model
Ai              =       0;          nAi             = size(Ai,1);
Bi1             =       [1 0];      mBil            = size(Bi1,2);
Bi2             =       -1;         mBi2            = size(Bi2,2);
Ci              =       [1;0];
Di1             =       [0 0 ; 0 1];    
Di2             =       [0;0];  

% Integrator state space system
sys_int_Az         =   ss(Ai,Bi1,Ci,Di1,'statename',{'eIAz'},'Inputname', {'Az_meas','q_meas'},'OutputName', {'eIAz_meas','q_meas'});

% Solve Dual ARE for Observer design model
[L, ~, ~] = lqr(Aw',Cw_meas',w,v);

% Define combined state space observer, controller and integrator systems
Aco              =       [Ai,zeros(nAi,nAw);L'*Ci,Aw - Bw*K_lqr - L'*Cw_meas];
Bco1             =       [Bi1;L'*Di1];
Bco2             =       [Bi2;-1;0;0];  mBco2    =   size(Bco2,2);
Cco              =       [0,-K_lqr];
Dco1             =       zeros(1,2);   nDco1    =   size(Dco1,1);
Dco2            =   0;
% Define compensator to determine loop gain at input
obsv_lqr_comp   =   ss(Aco,Bco1,Cco,Dco1,'stateName', {'eIAz_meas','eIAz_hat','alpha_hat','q_hat'},'InputName',{'Az_meas','q_meas'},'OutputName',{'dele'});

sys_mod = ss(Ap_sp,Bp_sp,[Zalpha 0; 0 1],[Zdele;0]);
sys_mod_act = sys_mod*sys_act;

    A_comp_cl = [Ap_sp Bp_sp zeros(2,5); zeros(2,2) Ap_act Bp_act*[0 -K_lqr];
                 Bco1*sys_mod_act.c  (Aco + Bco1*sys_mod_act.d*[0 -K_lqr])];

    B_comp_cl = [zeros(4,1);Bco2];
    C_comp_cl = [eye(4) zeros(4,4);zeros(1,6) 1 0;
                zeros(1,7) 1; Zalpha 0  Zdele zeros(1,5)];
    D_comp_cl = zeros(7,1);


eig_cl_obs_sys   =   squeeze(min(real(eig(A_comp_cl))));
min_closed_loop_eigs_sys_comp(jj) = eig_cl_obs_sys;
    sys_comp_cl = ss(A_comp_cl,B_comp_cl,C_comp_cl,D_comp_cl);
    [y_obs,~,x_obs_cl] = step(sys_comp_cl,t_sec);
    y_reg_obs = y_obs(:,7);
    
    % Percent undershoot
    percent_undershoot_obs = abs(min(y_reg_obs))*100;
    % Percent overshoot
    percent_overshoot_obs = ( abs(max(y_reg_obs))-1 )*100;    
    % Fin displacement
    max_elevator_angle_deg_obs = 180/pi*max(abs(x_obs_cl(:, 7)));
    % Fin displacement rate
    max_elevator_rate_dps_obs = 180/pi*max(abs(x_obs_cl(:, 8)));

    C_sen_obs = [Zalpha 0 Zdele zeros(1,5);
            0,1,0,0,0,0,0,0];

    sys_comp_cl_sen = ss(A_comp_cl,B_comp_cl,C_sen_obs,[0;0]);
    sys_comp_cl_sen_tf = tf(sys_comp_cl_sen);

    T_Az_obs = max(squeeze(mag2db(abs(freqresp(sys_comp_cl_sen_tf(1),w_rps))))); 
    tf_S_Az_obs = 1 - sys_comp_cl_sen_tf(1);
    S_Az_obs = max(squeeze(mag2db(abs(freqresp(tf_S_Az_obs,w_rps)))));

    T_q_obs = max(squeeze(mag2db(abs(freqresp(sys_comp_cl_sen_tf(2),w_rps)))));
    tf_S_q_obs = 1 - sys_comp_cl_sen_tf(2);
    S_q_obs = max(squeeze(mag2db(abs(freqresp(tf_S_q_obs,w_rps)))));


    compensator = ss(Aco, Bco1, Cco, Dco1);
    Lu_obs = compensator*sys_mod_act;[nLu_obs, mLu_obs] = size(Lu_obs.a);
    Lu_obs.c = -Lu_obs.c;
    Lu_obs.d= -Lu_obs.d;
    lu_tf_obs   =   tf(Lu_obs);
    % Evaluate loop gain at plant input over frequency range
    Lu_eval_obs = zeros(size(w_rps)); 
    eye_nLu_obs = eye(nLu_obs);
    for i = 1 :length(w_rps)
        resolvent = (1i*w_rps(i)*eye_nLu_obs - Lu_obs.a)\eye_nLu_obs;
        Lu_eval_obs(i) = Lu_obs.d + Lu_obs.c*resolvent*Lu_obs.b;
    end

    % Get stability margins
    [GM_obs, PM_obs, VM_obs] = stability_margins( w_rps, Lu_eval_obs, 'off');
    GM_obs_dB  = 20*log10(GM_obs);
    GM_obs_dB  = min(GM_obs_dB(GM_obs_dB > 0));
    PM_obs_deg = min(PM_obs(PM_obs > 0));

    % Compute crossover frequency of loop gain at plant input
    Lu_obs_dB = 20*log10(abs(Lu_eval_obs));
    wc_obs_rps = crosst(Lu_obs_dB, w_rps);
    wc_obs_Hz = wc_obs_rps/6.28;
    w_diff = (design_data(jj,7) - wc_obs_Hz);

    design_data_obs(jj,:) = [percent_overshoot_obs, percent_undershoot_obs,...
    max_elevator_angle_deg_obs, max_elevator_rate_dps_obs, wc_obs_Hz, GM_obs_dB, PM_obs_deg,...
    T_Az_obs, S_Az_obs, T_q_obs, S_q_obs, w_diff];
    fprintf('Compensator Design Iteration %2.0f complete.\n',jj);
end
sel_dd_obs_OS = 1;
sel_dd_obs_US = 2;
sel_dd_obs_findisp = 3;
sel_dd_obs_findisprate = 4;
sel_dd_obs_wc = 5;
sel_dd_obs_GM = 6;
sel_dd_obs_PM = 7;
sel_dd_obs_T_Az = 8;
sel_dd_obs_S_Az = 9;
sel_dd_obs_T_q = 10;
sel_dd_obs_S_q = 11;
sel_dd_obs_w_diff = 12;

figure(15);
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_OS),'*','linewidth',0.75); hold on
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_US),'r*','linewidth',0.75); hold on 
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel('[%]','fontsize', 14);
title('Design Chart: Percent Overshoot and Undershoot ','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 
h = legend('%OS','%US');
set(h,'location','best');

figure(16);
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_findisp)*644,'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel({'fin Displacement (Deg) ','[(20g ft/s^2 Az)]'},'fontsize', 14);
title('Design Chart: Fin displacement with Compensator','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(17);
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_findisprate)*644,'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel({'Fin Displacement Rate ( Deg/sec)','[(20g ft/s^2 Az)]'},'fontsize', 14);
title('Design Chart: Fin Displacement Rate with Compensator','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(18);
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_GM),'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel('Gain Margin [dB]','fontsize', 14);
title('Design Chart: GM at Plant Input with Observer','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(19);
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_PM),'*','linewidth',0.75); hold on
xlabel('Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel('Phase Margin [deg]','fontsize', 14);
title('Design Chart: Phase Margin at Plant Input with Compensator','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(20)
semilogx(q_obs(1,:), design_data_obs(:,sel_dd_obs_wc)*6.28,'ro','linewidth',0.75); hold on
xlabel('Penalty q_0', 'fontsize', 14);
ylabel('CrossOver Frequency','fontsize', 14);
title('Design Chart: Crossover frequency','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 
figure(21)
plot(design_data(:,sel_dd_wc), min_closed_loop_eigs_sysw,'*','linewidth',0.75); hold on
xlabel('Cross Over frequency (Hz)', 'fontsize', 14);
ylabel('Min(real(eig(Aw - Bw*Kc))','fontsize', 14);
title('Min Closed loop eigen values state feedback controller','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 
figure(22)
plot(design_data_obs(:,sel_dd_obs_wc), min_closed_loop_eigs_sys_comp,'*','linewidth',0.75); hold on
xlabel('Cross Over frequency (Hz)', 'fontsize', 14);
ylabel('Min closed loop eigen values','fontsize', 14);
title('Min Closed loop eigen values with compensator','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 

figure(23);
plot(design_data_obs(:,sel_dd_obs_wc), design_data_obs(:,sel_dd_obs_w_diff),'*',design_data_obs(23,sel_dd_obs_wc), design_data_obs(23,sel_dd_obs_w_diff),'ro','linewidth',0.75); hold on
xlabel('Compensator Crossover, \omega_c [Hz]', 'fontsize', 14);
ylabel({'\Delta_\omega','[(20g ft/s^2 Az)]'},'fontsize', 14);
title('Crossoverfrequecny difference','fontsize',14);
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on 


figure(24);
semilogx(w_rps1, design_data_obs(:, sel_dd_obs_S_Az),'b*','linewidth',0.75); hold on
semilogx(w_rps1, design_data_obs(:, sel_dd_obs_T_Az),'r*','linewidth',0.75); hold on
xlabel('Frequency, \omega rad/sec', 'fontsize', 14);
ylabel('| S | & | T | Az[dB]', 'fontsize', 14);
title('Sensitivity & Complimentary sensitivity of Az with observer','fontsize',14);
legend('|S|','|T|');
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on

figure(25);
semilogx(w_rps1, design_data(:, sel_dd_obs_S_q),'b*','linewidth',0.75); hold on
semilogx(w_rps1, design_data(:, sel_dd_obs_T_q),'r*','linewidth',0.75); hold on
xlabel('Frequency, \omega rad/sec', 'fontsize', 14);
ylabel('| S | & | T | q [dB]', 'fontsize', 14);
title('Sensitivity & Complimentary sensitivity of q with observer','fontsize',14);
legend('|S|','|T|');
set(gca,'fontsize', 14);
set(gcf,'color', 'w');
grid on

dp_select_obs = 19;
wc_design_obs = design_data_obs(dp_select_obs,sel_dd_obs_wc);
q_obs_dp = q_obs(dp_select_obs);
w_obs_dp = eye(3)*q_obs_dp;
v_obs_dp = v;
[L_dp, ~, ~] = lqr(Aw',Cw_meas',w,v);

% Define combined state space observer, controller and integrator systems
Aco_dp              =       [Ai,zeros(nAi,nAw);L_dp'*Ci,Aw - Bw*K_lqr - L_dp'*Cw_meas];
Bco1_dp             =       [Bi1;L_dp'*Di1];
Bco2_dp             =       [Bi2;-1;0;0];  mBco2_dp    =   size(Bco2,2);
Cco_dp              =       [0,-K_lqr];
Dco1_dp             =       zeros(1,2);   nDco1_dp    =   size(Dco1,1);
Dco2_dp             =   0;
% Define compensator to determine loop gain at input
obsv_lqr_comp_dp   =   ss(Aco_dp,Bco1_dp,Cco_dp,Dco1_dp,'stateName', {'eIAz_meas','eIAz_hat','alpha_hat','q_hat'},'InputName',{'Az_meas','q_meas'},'OutputName',{'dele'});



    A_comp_cl_dp = [Ap_sp Bp_sp zeros(2,5); zeros(2,2) Ap_act Bp_act*[0 -K_lqr];
                 Bco1*sys_mod_act.c  (Aco + Bco1*sys_mod_act.d*[0 -K_lqr])];

    B_comp_cl_dp = [zeros(4,1);Bco2];
    C_comp_cl_dp = [eye(4) zeros(4,4);zeros(1,6) 1 0;
                zeros(1,7) 1; Zalpha 0  Zdele zeros(1,5)];
    D_comp_cl_dp = zeros(7,1);

% Plant with actuator output data:
Cp_act_meas = Cp_sp_act([3,2],:);
Dp_act_meas = Dp_sp_act([3,2],:);
nAp_sp_act  = size(Ap_sp_act,1);
Z_dp        = eye(nDco1_dp) - Dco1_dp*Dp_act_meas;
Zinv_dp     = Z/eye(nDco1_dp);

compensator_dp = ss(Aco_dp, Bco1_dp, Cco_dp, Dco1_dp);
Lu_obs_dp = compensator_dp*sys_mod_act;[nLu_obs_dp, mLu_obs_dp] = size(Lu_obs_dp.a);
Lu_obs_dp.c = -Lu_obs_dp.c;
Lu_obs_dp.d= -Lu_obs_dp.d;

lu_tf_obs_dp   =   tf(Lu_obs_dp);
%% Bode diagram
figure(30)
margin(Lu_obs_dp)
set(gcf,'color','w');
C_comp_cl_dp_exe = [eye(8);C_sen_obs(1,:)];
D_comp_cl_dp_exe = zeros(9,1);
save('A_comp_cl_dp.MAT','A_comp_cl_dp');
save('B_comp_cl_dp.MAT','B_comp_cl_dp');
save('C_comp_cl_dp_exe.MAT','C_comp_cl_dp_exe');
save('D_comp_cl_dp_exe.MAT','D_comp_cl_dp_exe');