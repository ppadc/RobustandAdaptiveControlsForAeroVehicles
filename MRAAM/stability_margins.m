% FUNCTION stability_margins.m approximates SISO gain margin, phase margin, 
% and vector margin from frequency response data input pointwise with 
% corresponding frequency values. This function relies on basic formulas and 
% does not fully account for special exceptions that may lead to erroneous 
% results. The user should not solely rely on this script for robustness
% evaluation of closed loop feedback control systems.
%
% Input
%   w_rps     - frequency in radians per second [N dimension vector]
%   L         - frequency response data [N dimension vector]
%   plot_flag - 'on' or 'off' [optional input]
%
% Output
%   GM      - gain margins when the phase angle is -180 degrees in [dB]
%   PM      - phase margins when |L| = 1 in [deg]
%   VM      - minimum distance to the critical point over w_rps
%
% History: Written by B. Dickinson, June 2022.
%-------------------------------------------------------------------------------

function [GM, PM, VM] = stability_margins( w_rps, L, plot_flag)

nw = length(w_rps);

% Compute phase angle
phase_rad = atan2(imag(L), real(L));

% Compute magnitude
Lmag = sqrt( real(L).^2 + imag(L).^2 );

% Find frequency intervals where phase crosses -180 deg
phase_logic = phase_rad > 0;
w_phase_interval = find(phase_logic(2:nw) - phase_logic(1:nw-1));
n_phase_intervals = length(w_phase_interval);

% Find frequency intervals where the gain crosses 1
gain_logic = Lmag < 1;
w_gain_interval = find(gain_logic(2:nw) - gain_logic(1:nw-1));
n_gain_intervals = length(w_gain_interval);

% Approximate the phase margins with linear interpolation
if n_gain_intervals ~= 0
  PM = zeros(1,n_gain_intervals);
  for ii = 1 : n_gain_intervals
    
    % Get interval indices
    int_g_i   = w_gain_interval(ii);
    int_g_ip1 = w_gain_interval(ii)+1;
    
    % Interpolate to |L| = 1
    w_Lmag_1_rps = w_rps(int_g_i) + ...
      ( w_rps(int_g_ip1) - w_rps(int_g_i) )/ ...
      (  Lmag(int_g_ip1) -  Lmag(int_g_i) )* ...
      (                1 -  Lmag(int_g_i) );
    
    % Interpolate to phase at phase margin, evaluate phase margin
    theta_star = phase_rad(int_g_i) + ...
               ( phase_rad(int_g_ip1) - phase_rad(int_g_i) )/ ...
               (     w_rps(int_g_ip1) -     w_rps(int_g_i) )* ...
               (         w_Lmag_1_rps -     w_rps(int_g_i) );
    PM(ii) = (theta_star + pi)*180/pi;
      
  end
else
  PM = NA;
end

% Approximate gain margins with linear interpolation
% Approximate the phase margins with linear interpolation
if n_phase_intervals ~= 0
  GM = zeros(1,n_phase_intervals);
  for ii = 1 : n_phase_intervals
    
    % Get interval indices
    int_p_i   = w_phase_interval(ii);
    int_p_ip1 = w_phase_interval(ii)+1;
    
    % Interpolate to -180 degrees
    w_m180_rps = w_rps(int_p_i) + ...
               ( w_rps(int_p_ip1) -     w_rps(int_p_i) )/ ...
           ( phase_rad(int_p_ip1) - phase_rad(int_p_i) )* ...
           (                    0 - phase_rad(int_p_i) );
    
    % Interpolate to L at gain margin, evaluate gain margin
    Lmag_star = Lmag(int_p_i) + ...
              ( Lmag(int_p_ip1) -  Lmag(int_p_i) )/ ...
              (w_rps(int_p_ip1) - w_rps(int_p_i) )* ...
              (     w_m180_rps  - w_rps(int_p_i) );
    GM(ii) = 20*log10(1/Lmag_star);
      
  end
else
  GM = Inf;
end

% Evaluate return difference
RD = 1 + L;
RDmag = sqrt( real(RD).^2 + imag(RD).^2 );
VM = min(RDmag);

if strcmp(plot_flag, 'bode')

  figure
  semilogx( w_rps./6.28, 20*log10(Lmag), 'w', 'linewidth', 1);
  xlabel('Frequency [Hz]', 'fontsize', 16);
  ylabel('Magnitude [dB]', 'fontsize', 16);
  set(gca,'fontsize',16,'color','k','gridcolor','w','xcolor','y', ...
          'ycolor','y', 'position', [0.1300   0.2009   0.7750   0.7241]);
  set(gcf,'color','k', 'position', [1134    523    560    253]);
  grid on
  shg

  figure
  semilogx( w_rps./6.28, phase_rad*180/pi, 'w', 'linewidth', 1);
  xlabel('Frequency [Hz]', 'fontsize', 16);
  ylabel('Phase Angle [deg]', 'fontsize', 16);
  set(gca,'fontsize',16,'color','k','gridcolor','w','xcolor','y', ...
          'ycolor','y', 'position', [0.1300   0.2009   0.7750   0.7241]);
  set(gcf,'color','k', 'position', [1134    164    560    253]);
  grid on
  shg
  
end
