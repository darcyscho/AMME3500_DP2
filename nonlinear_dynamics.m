function ddq = nonlinear_dynamics(q, dq, tau)
% Inputs:
%   q   = [q1; q2; q3] joint positions (rad)
%   dq  = [dq1; dq2; dq3] joint velocities (rad/s)
%   tau = [tau1; tau2; tau3] joint torques (Nm)
% Output:
%   ddq = [ddq1; ddq2; ddq3] joint accelerations (rad/s^2)

% Unpack inputs
q1 = q(1); q2 = q(2); q3 = q(3);
dq1 = dq(1); dq2 = dq(2); dq3 = dq(3);
tau1 = tau(1); tau2 = tau(2); tau3 = tau(3);


% === ddq1 ===

ddq1 = (-1*dq1 - 0.225*-9.81*8.4*cos(q1) - 0.22^2*3.8*ddq2 - 0.22*-9.81*3.8*cos(q1 + q2) + 2.0*0.22*0.45*3.8*sin(q2)*dq1*dq2 + 0.22*0.45*3.8*sin(q2)*dq2^2 - 0.22*0.45*3.8*cos(q2)*ddq2 - 0.13^2*1.114*ddq2 - 0.13^2*1.114*ddq3 - 0.13*-9.81*1.114*cos(q1 + q2 + q3) + 2.0*0.13*0.45*1.114*sin(q2 + q3)*dq1*dq2 + 2.0*0.13*0.45*1.114*sin(q2 + q3)*dq1*dq3 + 0.13*0.45*1.114*sin(q2 + q3)*dq2^2 + 2.0*0.13*0.45*1.114*sin(q2 + q3)*dq2*dq3 + 0.13*0.45*1.114*sin(q2 + q3)*dq3^2 - 0.13*0.45*1.114*cos(q2 + q3)*ddq2 - 0.13*0.45*1.114*cos(q2 + q3)*ddq3 + 2.0*0.13*0.44*1.114*sin(q3)*dq1*dq3 + 2.0*0.13*0.44*1.114*sin(q3)*dq2*dq3 + 0.13*0.44*1.114*sin(q3)*dq3^2 - 2.0*0.13*0.44*1.114*cos(q3)*ddq2 - 0.13*0.44*1.114*cos(q3)*ddq3 - -9.81*0.45*3.8*cos(q1) - -9.81*0.45*1.114*cos(q1) - -9.81*0.44*1.114*cos(q1 + q2) + 2.0*0.45*0.44*1.114*sin(q2)*dq1*dq2 + 0.45*0.44*1.114*sin(q2)*dq2^2 - 0.45*0.44*1.114*cos(q2)*ddq2 - 0.44^2*1.114*ddq2 + 2.0*tau1)/(0.225^2*8.4 + 0.22^2*3.8 + 2.0*0.22*0.45*3.8*cos(q2) + 0.13^2*1.114 + 2.0*0.13*0.45*1.114*cos(q2 + q3) + 2.0*0.13*0.44*1.114*cos(q3) + 0.45^2*3.8 + 0.45^2*1.114 + 2.0*0.45*0.44*1.114*cos(q2) + 0.44^2*1.114);
ddq1_num = simplify(subs(ddq1, params, values));

% === ddq2 ===
ddq2 = (-1*dq2 - 0.22^2*3.8*ddq1 - 0.22*-9.81*3.8*cos(q1 + q2) - 0.22*0.45*3.8*sin(q2)*dq1^2 - 0.22*0.45*3.8*cos(q2)*ddq1 - 0.13^2*1.114*ddq1 - 0.13^2*1.114*ddq3 - 0.13*-9.81*1.114*cos(q1 + q2 + q3) - 0.13*0.45*1.114*sin(q2 + q3)*dq1^2 - 0.13*0.45*1.114*cos(q2 + q3)*ddq1 + 2.0*0.13*0.44*1.114*sin(q3)*dq1*dq3 + 2.0*0.13*0.44*1.114*sin(q3)*dq2*dq3 + 0.13*0.44*1.114*sin(q3)*dq3^2 - 2.0*0.13*0.44*1.114*cos(q3)*ddq1 - 0.13*0.44*1.114*cos(q3)*ddq3 - -9.81*0.44*1.114*cos(q1 + q2) - 0.45*0.44*1.114*sin(q2)*dq1^2 - 0.45*0.44*1.114*cos(q2)*ddq1 - 0.44^2*1.114*ddq1 + 2.0*tau2)/(0.22^2*3.8 + 0.13^2*1.114 + 2.0*0.13*0.44*1.114*cos(q3) + 0.44^2*1.114);
ddq2_num = simplify(subs(ddq2, params, values));

% === ddq3 ===
ddq3 = (-1*dq3 - 0.13^2*1.114*(ddq1 + ddq2) - 0.13*1.114*(-9.81*cos(q1 + q2 + q3) + 0.45*sin(q2 + q3)*dq1^2 + 0.45*cos(q2 + q3)*ddq1 + 0.44*sin(q3)*dq1^2 + 2.0*0.44*sin(q3)*dq1*dq2 + 0.44*sin(q3)*dq2^2 + 0.44*cos(q3)*ddq1 + 0.44*cos(q3)*ddq2) + 2.0*tau3)/(0.13^2*1.114);
ddq3_num = simplify(subs(ddq3, params, values));

% Output
ddq = [ddq1_num; ddq2_num; ddq3_num];
end
