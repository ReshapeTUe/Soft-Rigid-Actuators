clear all; close all;clc;

%% Set text interperters to 'latex'
set(0, 'defaultTextInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

%% Plot Settings
fsize = 15; % Font size
lw = 3;   % Line width

%% Global Parameters
global m d_th d g l0 K k e_N

Tspan=[0 20];
m = 0.2;
d_th = 0.1;

g = 9.81;
l0 = 0.10;
K = 30;

e_N = 0;


q_init_r= [0; 0];% Initial condition (should satisfy the constraints)
q_init_f= [l0; 0; 0; 0];% Initial condition (should satisfy the constraints)

options_r = odeset('Events',@(t,q) ImpactRigid(t,q), 'RelTol', 1e-6, 'AbsTol', 1e-9); % For detecting impacts in ODE solver
options_f = odeset('Events',@(t,q) ImpactFlex(t,q), 'RelTol', 1e-6, 'AbsTol', 1e-9);  % For detecting impacts in ODE solver

%% Model Rigid
NDOF = 1;
[tr,qr, ~, ~, Impact] = ode45(@(t,q) SIMrigid(t,q),Tspan,q_init_r,options_r);

if Impact
    while tr(end) < Tspan(end)      % Keeps loop running for as long as events keep interupping the simulation
        qminus = qr(end,1:NDOF)';           % Ante-impact joint positions
        qdminus = qr(end,NDOF+1:2*NDOF)';   % Ante-impact joint velocities
        J_c_N_imp = cos(qminus)*l0;
        M = (1/3)*m*l0^2;
        ImpactMat = [M -J_c_N_imp'; J_c_N_imp 0];  % Invertible matrix structured as [M, J_N^T ; J_N, 0]
        ImpactVec = [M*qdminus; e_N*J_c_N_imp*qdminus];               % Vector Structured as [M*q^-; e_N*J_N*q^-]
        qdplusImpulse = ImpactMat\ImpactVec;   % Calculates post-impact joint velocities and Normal Impulses from contact
        qplus = qminus;                 % Joint positions, motor positions and motor velocities remain constant
        qplus(2) = qdplusImpulse(1);  % Update old joint velocities with new joint velocities
        
        % Run Simulation After Event and Merge Data
        [tr_t,qr_t] = ode45(@(t,q) SIMrigid(t,q),[tr(end) Tspan(end)],qplus,options_r);
        tr = [tr; tr_t];
        qr = [qr; qr_t];
    end
end

F_r = zeros(length(tr),1);

for i = 1:length(tr)
    U = -tr(i)/20*0.5*pi;
    q      = qr(i,1:NDOF)';        % Joint positions  [rad]
    qd     = qr(i,NDOF+1:2*NDOF)'; % Joint velocities [rad/s]

    M = (1/3)*m*l0^2;
    C = 0.5*m*g*l0*cos(q)+K*(q-U)+d_th*qd;
    J_c_N_imp = cos(q)*l0;
    Jd_c_N_imp = -sin(q)*l0*qd;

    if m*g*sin(q)*l0-m*g*sin(-0.25*pi)*l0 <=0  && qd*cos(q)*l0<=0
        A = [M, -J_c_N_imp';
            J_c_N_imp,    0];
        b = [-C;
            -Jd_c_N_imp*qd];
        z = linsolve(A,b); % Joint accelerations and contact forces
    else
        A = M;
        b = -C;
        z=[0; 0];
    end

    g_r(i) = sin(q)*l0 - sin(-0.25*pi)*l0; 

    F_r(i) = z(end);

end

%% Model Flexible (Least Stiff)

NDOF = 2;

k = 500;
d = k/1000;


[tf,qf, ~, ~, Impact] = ode45(@(t,q) SIMflex(t,q),Tspan,q_init_f,options_f);

if Impact
    while tf(end) < Tspan(end)
        qminus = qf(end,1:NDOF)';
        qdminus = qf(end,NDOF+1:2*NDOF)';   % Ante-impact joint velocities
        J_c_N_imp = [sin(qminus(2)) cos(qminus(2))*qminus(1)];
        M = [0.25*m 0;
            0 (1/3)*m*qminus(1)^2];
        ImpactMat = [M -J_c_N_imp'; J_c_N_imp 0];  % Invertible matrix structured as [M, J_N^T ; J_N, 0]
        ImpactVec = [M*qdminus; e_N*J_c_N_imp*qdminus];               % Vector Structured as [M*q^-; e_N*J_N*q^-]
        qdplusImpulse = ImpactMat\ImpactVec;   % Calculates post-impact joint velocities and Normal Impulses from contact
        qplus = qminus;                % Joint positions, motor positions and motor velocities remain constant
        qplus(NDOF+1:2*NDOF) = qdplusImpulse(1:NDOF);  % Update old joint velocities with new joint velocities

        [tf_t,qf_t] = ode45(@(t,q) SIMflex(t,q),[tf(end) Tspan(end)],qplus,options_f);
 
        tf = [tf;tf_t];
        qf = [qf;qf_t];
    end
end

F_f = zeros(length(tf),1);

for i = 1:length(tf)
    U = -tf(i)/20*0.5*pi;
    q      = qf(i,1:NDOF)';        % Joint positions  [rad]
    qd     = qf(i,NDOF+1:2*NDOF)'; % Joint velocities [rad/s]

    M = [0.25*m 0;
     0 (1/3)*m*q(1)^2];

    C = [-(1/3)*m*q(1)*qd(2)^2 + 0.5*m*g*sin(q(2)) + k*(q(1)-l0) + d*qd(1);
    (2/3)*m*qd(1)*q(1)*qd(2)+ 0.5*m*g*cos(q(2))*q(1)+K*(q(2)-U)+d_th*qd(2)];

    J_c_N_imp = [sin(q(2)) cos(q(2))*q(1)];
    Jd_c_N_imp = [cos(q(2))*qd(2)   -sin(q(2))*q(1)*qd(2)+cos(q(2))*qd(1) ];


    if m*g*sin(q(2))*q(1)-m*g*sin(-0.25*pi)*l0 <= 0 && qd(2)*cos(q(2))*q(1) + sin(q(2))*qd(1) <= 0
        A = [M, -J_c_N_imp';
            J_c_N_imp,    0];
        b = [-C;
            -Jd_c_N_imp*qd];
        z = linsolve(A,b); % Joint accelerations and contact forces
    else
        A = M;
        b = -C;
        z = [0;0;0];
    end
    
    g_f(i) = sin(q(2))*q(1) - sin(-0.25*pi)*l0;

    F_f(i) = z(end);
end
g_1 = g_f;
q_1 = qf;
t_1 = tf;
F_1 = F_f;


%% Model Flexible (mid Stiff)

NDOF = 2;

k = 5000;
d = k/1000;


[tf,qf, ~, ~, Impact] = ode45(@(t,q) SIMflex(t,q),Tspan,q_init_f,options_f);

if Impact
    while tf(end) < Tspan(end)
        qminus = qf(end,1:NDOF)';
        qdminus = qf(end,NDOF+1:2*NDOF)';   % Ante-impact joint velocities
        J_c_N_imp = [sin(qminus(2)) cos(qminus(2))*qminus(1)];
        M = [0.25*m 0;
            0 (1/3)*m*qminus(1)^2];
        ImpactMat = [M -J_c_N_imp'; J_c_N_imp 0];  % Invertible matrix structured as [M, J_N^T ; J_N, 0]
        ImpactVec = [M*qdminus; e_N*J_c_N_imp*qdminus];               % Vector Structured as [M*q^-; e_N*J_N*q^-]
        qdplusImpulse = ImpactMat\ImpactVec;   % Calculates post-impact joint velocities and Normal Impulses from contact
        qplus = qminus;                % Joint positions, motor positions and motor velocities remain constant
        qplus(NDOF+1:2*NDOF) = qdplusImpulse(1:NDOF);  % Update old joint velocities with new joint velocities

        [tf_t,qf_t] = ode45(@(t,q) SIMflex(t,q),[tf(end) Tspan(end)],qplus,options_f);
 
        tf = [tf;tf_t];
        qf = [qf;qf_t];
    end
end

F_f = zeros(length(tf),1);

for i = 1:length(tf)
    U = -tf(i)/20*0.5*pi;
    q      = qf(i,1:NDOF)';        % Joint positions  [rad]
    qd     = qf(i,NDOF+1:2*NDOF)'; % Joint velocities [rad/s]

    M = [0.25*m 0;
     0 (1/3)*m*q(1)^2];

    C = [-(1/3)*m*q(1)*qd(2)^2 + 0.5*m*g*sin(q(2)) + k*(q(1)-l0) + d*qd(1);
    (2/3)*m*qd(1)*q(1)*qd(2)+ 0.5*m*g*cos(q(2))*q(1)+K*(q(2)-U)+d_th*qd(2)];

    J_c_N_imp = [sin(q(2)) cos(q(2))*q(1)];
    Jd_c_N_imp = [cos(q(2))*qd(2)   -sin(q(2))*q(1)*qd(2)+cos(q(2))*qd(1) ];


    if m*g*sin(q(2))*q(1)-m*g*sin(-0.25*pi)*l0 <= 0 && qd(2)*cos(q(2))*q(1) + sin(q(2))*qd(1) <= 0
        A = [M, -J_c_N_imp';
            J_c_N_imp,    0];
        b = [-C;
            -Jd_c_N_imp*qd];
        z = linsolve(A,b); % Joint accelerations and contact forces
    else
        A = M;
        b = -C;
        z = [0;0;0];
    end
    
    g_f(i) = sin(q(2))*q(1) - sin(-0.25*pi)*l0;

    F_f(i) = z(end);
end
g_2 = g_f;
q_2 = qf;
t_2 = tf;
F_2 = F_f;


%% Model Flexible (High Stiff)

NDOF = 2;

k = 50000;
d = k/1000;


[tf,qf, ~, ~, Impact] = ode45(@(t,q) SIMflex(t,q),Tspan,q_init_f,options_f);

if Impact
    while tf(end) < Tspan(end)
        qminus = qf(end,1:NDOF)';
        qdminus = qf(end,NDOF+1:2*NDOF)';   % Ante-impact joint velocities
        J_c_N_imp = [sin(qminus(2)) cos(qminus(2))*qminus(1)];
        M = [0.25*m 0;
            0 (1/3)*m*qminus(1)^2];
        ImpactMat = [M -J_c_N_imp'; J_c_N_imp 0];  % Invertible matrix structured as [M, J_N^T ; J_N, 0]
        ImpactVec = [M*qdminus; e_N*J_c_N_imp*qdminus];               % Vector Structured as [M*q^-; e_N*J_N*q^-]
        qdplusImpulse = ImpactMat\ImpactVec;   % Calculates post-impact joint velocities and Normal Impulses from contact
        qplus = qminus;                % Joint positions, motor positions and motor velocities remain constant
        qplus(NDOF+1:2*NDOF) = qdplusImpulse(1:NDOF);  % Update old joint velocities with new joint velocities

        [tf_t,qf_t] = ode45(@(t,q) SIMflex(t,q),[tf(end) Tspan(end)],qplus,options_f);
 
        tf = [tf;tf_t];
        qf = [qf;qf_t];
    end
end

F_f = zeros(length(tf),1);

for i = 1:length(tf)
    U = -tf(i)/20*0.5*pi;
    q      = qf(i,1:NDOF)';        % Joint positions  [rad]
    qd     = qf(i,NDOF+1:2*NDOF)'; % Joint velocities [rad/s]

    M = [0.25*m 0;
     0 (1/3)*m*q(1)^2];

    C = [-(1/3)*m*q(1)*qd(2)^2 + 0.5*m*g*sin(q(2)) + k*(q(1)-l0) + d*qd(1);
    (2/3)*m*qd(1)*q(1)*qd(2)+ 0.5*m*g*cos(q(2))*q(1)+K*(q(2)-U)+d_th*qd(2)];

    J_c_N_imp = [sin(q(2)) cos(q(2))*q(1)];
    Jd_c_N_imp = [cos(q(2))*qd(2)   -sin(q(2))*q(1)*qd(2)+cos(q(2))*qd(1) ];


    if m*g*sin(q(2))*q(1)-m*g*sin(-0.25*pi)*l0 <= 0 && qd(2)*cos(q(2))*q(1) + sin(q(2))*qd(1) <= 0
        A = [M, -J_c_N_imp';
            J_c_N_imp,    0];
        b = [-C;
            -Jd_c_N_imp*qd];
        z = linsolve(A,b); % Joint accelerations and contact forces
    else
        A = M;
        b = -C;
        z = [0;0;0];
    end
    
    g_f(i) = sin(q(2))*q(1) - sin(-0.25*pi)*l0;

    F_f(i) = z(end);
end
q_3 = qf;
t_3 = tf;
F_3 = F_f;
g_3 = g_f;
%% Model Flexible (Highest Stiff)

NDOF = 2;

k = 5000000;
d = k/1000;


[tf,qf, ~, ~, Impact] = ode45(@(t,q) SIMflex(t,q),Tspan,q_init_f,options_f);

if Impact
    while tf(end) < Tspan(end)
        qminus = qf(end,1:NDOF)';
        qdminus = qf(end,NDOF+1:2*NDOF)';   % Ante-impact joint velocities
        J_c_N_imp = [sin(qminus(2)) cos(qminus(2))*qminus(1)];
        M = [0.25*m 0;
            0 (1/3)*m*qminus(1)^2];
        ImpactMat = [M -J_c_N_imp'; J_c_N_imp 0];  % Invertible matrix structured as [M, J_N^T ; J_N, 0]
        ImpactVec = [M*qdminus; e_N*J_c_N_imp*qdminus];               % Vector Structured as [M*q^-; e_N*J_N*q^-]
        qdplusImpulse = ImpactMat\ImpactVec;   % Calculates post-impact joint velocities and Normal Impulses from contact
        qplus = qminus;                % Joint positions, motor positions and motor velocities remain constant
        qplus(NDOF+1:2*NDOF) = qdplusImpulse(1:NDOF);  % Update old joint velocities with new joint velocities

        [tf_t,qf_t] = ode45(@(t,q) SIMflex(t,q),[tf(end) Tspan(end)],qplus,options_f);
 
        tf = [tf;tf_t];
        qf = [qf;qf_t];
    end
end

F_f = zeros(length(tf),1);

for i = 1:length(tf)
    U = -tf(i)/20*0.5*pi;
    q      = qf(i,1:NDOF)';        % Joint positions  [rad]
    qd     = qf(i,NDOF+1:2*NDOF)'; % Joint velocities [rad/s]

    M = [0.25*m 0;
     0 (1/3)*m*q(1)^2];

    C = [-(1/3)*m*q(1)*qd(2)^2 + 0.5*m*g*sin(q(2)) + k*(q(1)-l0) + d*qd(1);
    (2/3)*m*qd(1)*q(1)*qd(2)+ 0.5*m*g*cos(q(2))*q(1)+K*(q(2)-U)+d_th*qd(2)];


    J_c_N_imp = [sin(q(2)) cos(q(2))*q(1)];
    Jd_c_N_imp = [cos(q(2))*qd(2)   -sin(q(2))*q(1)*qd(2)+cos(q(2))*qd(1) ];


    if m*g*sin(q(2))*q(1)-m*g*sin(-0.25*pi)*l0 <= 0 && qd(2)*cos(q(2))*q(1) + sin(q(2))*qd(1) <= 0
        A = [M, -J_c_N_imp';
            J_c_N_imp,    0];
        b = [-C;
            -Jd_c_N_imp*qd];
        z = linsolve(A,b); % Joint accelerations and contact forces
    else
        A = M;
        b = -C;
        z = [0;0;0];
    end
    
    g_f(i) = sin(q(2))*q(1) - sin(-0.25*pi)*l0;

    F_f(i) = z(end);
end



q_4 = qf;
t_4 = tf;
F_4 = F_f;
g_4 = g_f;

%% Plots

t_offset = 0;

Fig1 = figure;
box on, grid on, hold on
plot(tr-t_offset,F_r,'LineWidth',lw, 'Color',[0, 0, 128]/(255))
plot(t_4-t_offset,F_4,'--','LineWidth',lw, 'Color',[255, 165, 0]/(255))
plot(t_3-t_offset,F_3,'--','LineWidth',lw, 'Color',[69, 139, 0]/(255))
plot(t_2-t_offset,F_2,'--','LineWidth',lw, 'Color',[205, 38, 38]/(255))
plot(t_1-t_offset,F_1,'--','LineWidth',lw, 'Color',[47, 79, 79]/(255))
xline(10-t_offset, 'k:','LineWidth',lw);
xlim([0 Tspan(end)-t_offset])
ax = gca;
ax.FontSize = fsize; 
ylabel('Contact Force $\lambda_N$ (N)', 'FontSize', fsize)
xlabel('Time $t$ (s)', 'FontSize', fsize)
ylim([0 400])
legend('Rigid-body ($k_{x} = \infty$)', 'Flexible-body ($k_{x} = k_{x,max}$)', 'Flexible-body ($k_{x} = 0.01k_{x,max}$)' ,'Flexible-body ($k_{x} = 0.001 k_{x,max}$)','Flexible-body ($k_{x} = 0.0001 k_{x,max}$)' , 'FontSize', fsize, 'Location', 'northwest')



Fig2 = figure;
box on, grid on, hold on
plot(tr-t_offset,qr(:,1),'LineWidth',lw, 'Color',[0, 0, 128]/(255))
plot(t_4-t_offset,q_4(:,2),'--','LineWidth',lw, 'Color',[255, 165, 0]/(255))
plot(t_3-t_offset,q_3(:,2),'--','LineWidth',lw, 'Color',[69, 139, 0]/(255))
plot(t_2-t_offset,q_2(:,2),'--','LineWidth',lw, 'Color',[205, 38, 38]/(255))
plot(t_1-t_offset,q_1(:,2),'--','LineWidth',lw, 'Color',[47, 79, 79]/(255))
xline(10-t_offset, 'k:','LineWidth',lw);
xlim([0 Tspan(end)-t_offset])
ax = gca;
ax.FontSize = fsize; 
ylabel('Joint Rotation $\theta$ (rad)', 'FontSize', fsize)
xlabel('Time $t$ (s)', 'FontSize', fsize)


Fig3 = figure;
box on, grid on, hold on
plot(Tspan,[l0 l0],'LineWidth',lw, 'Color',[0, 0, 128]/(255))
plot(t_4-t_offset,q_4(:,1),'--','LineWidth',lw, 'Color',[255, 165, 0]/(255))
plot(t_3-t_offset,q_3(:,1),'--','LineWidth',lw, 'Color',[69, 139, 0]/(255))
plot(t_2-t_offset,q_2(:,1),'--','LineWidth',lw, 'Color',[205, 38, 38]/(255))
plot(t_1-t_offset,q_1(:,1),'--','LineWidth',lw, 'Color',[47, 79, 79]/(255))
xline(10-t_offset, 'k:','LineWidth',lw);
xlim([0 Tspan(end)-t_offset])
ax = gca;
ax.FontSize = fsize; 
ylabel('Beam Length $x$ (m)', 'FontSize', fsize)
xlabel('Time $t$ (s)', 'FontSize', fsize)



Fig4 = figure;
box on, grid on, hold on
plot(tr-t_offset,g_r,'LineWidth',lw, 'Color',[0, 0, 128]/(255))
plot(t_4-t_offset,g_4,'--','LineWidth',lw, 'Color',[255, 165, 0]/(255))
plot(t_3-t_offset,g_3,'--','LineWidth',lw, 'Color',[69, 139, 0]/(255))
plot(t_2-t_offset,g_2,'--','LineWidth',lw, 'Color',[205, 38, 38]/(255))
plot(t_1-t_offset,g_1,'--','LineWidth',lw, 'Color',[47, 79, 79]/(255))
xline(10-t_offset, 'k:','LineWidth',lw);
xlim([0 Tspan(end)-t_offset])
ax = gca;
ax.FontSize = fsize; 
ylabel('Distance to Surface $x$ (m)', 'FontSize', fsize)
xlabel('Time $t$ (s)', 'FontSize', fsize)


