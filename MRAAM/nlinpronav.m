

function dy = nlinpronav(t, y, He_rad, Np,A,B,C,D, nT)

% Define pointers to state variables
%--------------------------------------------------------------------------
% Pointers to states
sel_beta = 1;
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;

sel_alpha = 10;
sel_q = 11;
sel_u = 12;
sel_udot = 13;
sel_ei = 14;
sel_ei_hat = 15;
sel_alpha_hat = 16;
sel_qhat = 17;

% Preallocate solution vector
%--------------------------------------------------------------------------
dy = zeros(17,1);
% Preliminary terms to compute DE RHS
%--------------------------------------------------------------------------

% target and missile velocity magnitudes
VT = sqrt( y(sel_VT1)^2 + y(sel_VT2)^2 );

% relative positions and velocities
RTM1 = y(sel_RT1) - y(sel_RM1);
RTM2 = y(sel_RT2) - y(sel_RM2);
VTM1 = y(sel_VT1) - y(sel_VM1);
VTM2 = y(sel_VT2) - y(sel_VM2);

% relative distance
RTM = sqrt(RTM1^2 + RTM2^2);

% line of sight angle and time derivative
lambda     = atan2( RTM2, RTM1 );
lambda_dot = (RTM1*VTM2 - RTM2*VTM1)/RTM^2;

% closing velocity
VC = -(RTM1*VTM1 + RTM2*VTM2)/RTM;

% command acceleration (True Proportional Navigation)
nc = Np*VC*lambda_dot;

states = [y(sel_alpha); y(sel_q); y(sel_u); y(sel_udot); y(sel_ei); y(sel_ei_hat); y(sel_alpha_hat); y(sel_qhat)];

states_dot =  A*states + B*nc;
% DE RHS computations y = [beta, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]
%--------------------------------------------------------------------------
dy(1) =  nT/VT;
dy(2) = -VT*cos(y(sel_beta));
dy(3) =  VT*sin(y(sel_beta));
dy(4) =  y(sel_VM1);
dy(5) =  y(sel_VM2);
dy(6) =  nT*sin(y(sel_beta));
dy(7) =  nT*cos(y(sel_beta));
dy(8) = -nc*sin(lambda);
dy(9) =  nc*cos(lambda);

dy(10) = states_dot(1);
dy(11) = states_dot(2);
dy(12) = states_dot(3);
dy(13) = states_dot(4);
dy(14) = states_dot(5);
dy(15) = states_dot(6);
dy(16) = states_dot(7);
dy(17) = states_dot(8);

end