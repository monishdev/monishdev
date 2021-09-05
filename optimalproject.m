%Optimal Control Project - EEE 587
%Monish Dev Sudhakhar & Sai Praveen Daddala (TEAM 5)
%OPTIMAL CRUISE CONRTOL
clear all
clc
%Set up time vector
t = 0:0.1:0.45;
dt = 0.1;
%t_dis = 0:dt:t;

%Initialize SS Matrices
A = [0 -1 0 1; 0 0 0 0; 0 0 0 -1; 0 0 0.9 -1.2];
B = [0; 1; 0; 0];
%C = [1 1; 0 0]
Ak=[0 -1 0 1; 0 0 0 0; 0 0 0 -1; 0 0 0.6 -1.2];
B1 = [0; 1; 0; 0];
t1 = 0.4
Ak = exp(Ak*t1)

Bk1 = (exp(Ak*0.35)*B1)
Bk2 = (exp(Ak*0.15)*B1)
%UPDATED STATE VAR

E = [Ak Bk2; 0 0 0 0 0]
P = [Bk1;1]
U = eye(5)
O = [0; 0; 0; 0; 0]

% %Discretize SS Matrices
sysc = ss(E,P,U,O);
sysd = c2d(sysc, dt);

%A and B matrices
E = sysd.a;
P = sysd.b;

%A = eye(2) + A*dt;
%B = B*dt;

%Initialize LQR gains
Q = [1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
R = 1;
Q1 = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0 ];

Q1 = Q1*dt
R = R*dt

%Check for controllability
control = ctrb(E,P);
if rank(control) == 4
    disp("System is controllable")
else
    disp("System is not controllable")
end

%Set N = length of t
N = length(t);
S = Q1


%Initialize vectors
L = zeros(N-1,5);


%Solve recurrence equations for F
for K = N-1:-1:1
    L(K,:) = (inv(P' * S * P + R))*P' * S * E;
    S = E'* S * E + Q1 - L(K)' * P' * S * E;
    V =  E' * S * E
    H =  E' * S * E + Q1
    velocity = 1/t1 *(H)
end
for K = N-2:-1:1
    L(K,:) = (inv(P' * S * P + R))*P' * S * E;
    S = E'* S * E + Q1 - L(K)' * P' * S * E;
    V =  E' * S * E
    H1 =  E' * S * E + Q1
    velocity1 =1/t1 * (H1)
end
%Plotting headway & velocity
figure(1)

plot(t(1,:),H(:,5),'r--s')
hold on
plot(t(1,:),H1(:,5),'b:o')
legend('Headway for N-1','Headway for N-2')
title('Headway Profile')
xlabel('Time')
ylabel('Headway')
figure(2)
plot(t(1,:),velocity(:,5),'r--s')
hold on
plot(t(1,:),velocity1(:,5),'b:o')
ylim([0 5])
legend('Velocity for N-1','Velocity for N-2')
title('Velocity Profile')
xlabel('Time')
ylabel('Velocity')


%X(:,1) = [1;0];
%X_dlqr(:,1) = [1;0];
%P_dlqr = dare(A,B,Q,R);
%K_dlqr = inv(R)*B'*P_dlqr;



%Set initial conditions
X0 = [1 2 3 4 5];
X = zeros(5, N-2);
X(:,1) = X0;

%Use numerical integration to plot x
for K = 1:N-1
    u(K+1) = -L(K,:)*X(:,K);
    X(:,K+1) = E*X(:,K) + P*L(K,:)*X(:,K);
    %X(:,K) = X_dlqr*dt + X(:,K-1);
end

for K = 2:N-2
    u1(K+2) = -L(K,:)*X(:,K)
    X1(:,K+2) = E*X(:,K) + P*L(K,:)*X(:,K);
    %X(:,K) = X_dlqr*dt + X(:,K-1);
end
% Plotting Control and Position
figure(3)
plot(t,u,'m-^','LineWidth',1)
hold on
plot(t,u1,'r--d','LineWidth',2)
grid on
legend('Control profile for N-1','Control Profile for N-2')
title('Control Profile')
xlabel('time')
ylabel('Control Signal')
disp(X)
disp(X1)
figure(4)
plot(t,X(:,5),'r--s','LineWidth',2)
hold on
plot(t,X1(:,5),'g--o','LineWidth',2)
grid on
legend('Position Vector for k = N-1','Position Vector for k = N-2')
xlabel('time')
ylabel('Position')