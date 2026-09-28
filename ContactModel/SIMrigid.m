function xdot=SIMrigid(t,x)
global m d_th g l0 K


U = -t/20*0.5*pi;

NDOF = 1;

q      = x(1:NDOF);        % Joint positions  [rad]
qd     = x(NDOF+1:2*NDOF); % Joint velocities [rad/s]

M = (1/3)*m*l0^2;
Cg = 0.5*m*g*l0*cos(q)+K*(q-U)+d_th*qd;
J_c_N_imp = cos(q)*l0;
Jd_c_N_imp = -sin(q)*l0*qd;


if  sin(q)*l0-sin(-0.25*pi)*l0 <=0 && qd*cos(q)*l0<=0
    A = [M, -J_c_N_imp';
        J_c_N_imp,    0];
    b = [-Cg;
        -Jd_c_N_imp*qd]; 
else
    A = M;
    b = -Cg;
end

z = linsolve(A,b); % Joint accelerations and contact forces

% Compute time derivative of states
xdot = [qd;      
    z(1:NDOF)]; 

end



