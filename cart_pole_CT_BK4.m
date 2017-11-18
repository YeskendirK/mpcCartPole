%clear;
%clc;

%Constants
m=1; l=1; M=1; g=9.8;
A=[0 0 1 0;0 0 0 1;(m+M)*g/l/m 0 0 0;-m*g/M 0 0 0];
B=[0; 0; -1/l/M; 1/M];

%Objectives and constraints
IC=[pi/6;1;0;0];
zmax=10; zmin=-10;
umax=100; umin=-100;
thetamax=pi/2; thetamin=-pi/2;
dt=0.05;
h = 0.01;
t=0:h:10;
Q = 1*eye(4); R = 1*eye(1);
Qhalf = sqrtm(Q); Rhalf = sqrtm(R);

%Model Predictive Control
N=50; x=zeros(4,size(t,2)); u=zeros(1,size(t,2));
init=IC; umpc =0; x(:,1)=init;
for k=1:size(t,2)-1
    if rem(k,5) == 1
        cvx_begin 
            variables X(4,N+1) U(1,N)
            max(X(2,:)') <= zmax; min(X(2,:)') >= zmin;
            max(X(1,:)') <= thetamax; min(X(1,:)') >= thetamin;
            max(U') <= umax; min(U') >= umin;
            X(:,2:N+1) == dt*(A*X(:,1:N)+B*U)+X(:,1:N);
            X(:,1) == init; 
            X(:,N+1) == 0;
            minimize (norm([Qhalf*X(:,1:N); Rhalf*U],'fro'));
        cvx_end
        u(1,k)=U(1,1);
        umpc =U(1,1);
    end;
        k1 = dx(t(k+1),x(:,k),umpc);
        k2 = dx(t(k+1)+h/2, x(:,k)+0.5*h*k1,umpc);
        k3 = dx(t(k+1)+h/2, x(:,k)+0.5*h*k2,umpc);
        k4 = dx(t(k+1)+h,x(:,k)+h*k3,umpc);
        x(:,k+1) = x(:,k)+(k1+2*k2+2*k3+k4)*h/6;
        init = x(:,k+1);
end


save CartPoleReults.mat


figure;
subplot(2,1,1);
set(gca,'Fontsize',16);
stairs(t,[x(1,1:size(t,2))],'k'); 
ylabel('theta');
title('mpc');
subplot(2,1,2); 
set(gca,'Fontsize',16);
stairs(t,[x(2,1:size(t,2))],'k'); 
xlabel('t'); ylabel('x');

figure;
subplot(2,1,1);
set(gca,'Fontsize',16);
stairs(t,[x(3,1:size(t,2))],'k');2]); 
ylabel('theta_dot');
title('mpc');
subplot(2,1,2); 
set(gca,'Fontsize',16);
stairs(t,[x(4,1:size(t,2))],'k'); 
xlabel('t'); ylabel('x_dot');

figure;
stairs(t,[u(1,1:size(t,2))],'k');
xlabel('t'); ylabel('u');