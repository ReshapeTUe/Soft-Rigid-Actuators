function xdot=SIMflex(t,x)
global m d_th d g l0 K k


U = -t/20*0.5*pi;

NDOF = 2;

q      = x(1:NDOF);        % Joint positions  [rad]
qd     = x(NDOF+1:2*NDOF); % Joint velocities [rad/s]

M = [0.25*m 0;
     0 (1/3)*m*q(1)^2];

Cg = [-(1/3)*m*q(1)*qd(2)^2 + 0.5*m*g*sin(q(2)) + k*(q(1)-l0) + d*qd(1);
    (2/3)*m*qd(1)*q(1)*qd(2)+ 0.5*m*g*cos(q(2))*q(1)+K*(q(2)-U)+d_th*qd(2)];

J_c_N_imp = [sin(q(2)) cos(q(2))*q(1)];
Jd_c_N_imp = [cos(q(2))*qd(2)   -sin(q(2))*q(1)*qd(2)+cos(q(2))*qd(1) ];

if sin(q(2))*q(1)-sin(-0.25*pi)*l0 <= 0   && qd(2)*cos(q(2))*q(1) + sin(q(2))*qd(1) <= 0
    A = [M, -J_c_N_imp';
        J_c_N_imp,    0];
    b = [-Cg;
        -Jd_c_N_imp*qd]; 
    
else
    A = M;
    b = -Cg;
end

z = linsolve(A,b); % Joint accelerations and contact forces


  xdot = [qd;      
    z(1:NDOF)]; 


end

