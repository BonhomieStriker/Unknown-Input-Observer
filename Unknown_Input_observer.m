clear all
close all
clc
%% Iteration parameters
iter = 200;
n = 1:iter;
t = n * 0.01;

%% system parameters
A = [0.5 1
    0 -0.5];
B = [1 0
    0 1];
E = [1 1/2]';
C = [1 -1
    0 1];
L = [1 1/2
    0 1/4];    

%% Unknown Input Observer parameters
H = E / (C * E); %E - HCE = 0
T = (eye(2) - H * C) * B;%B - T - HCB = 0;
K1 =[-1 -2
    -0.5 -3];
N = A - H * C * A - K1 * C;%A - HCA - K1C - N = 0
K2 = N * H; 
K = K1 + K2;

%% Initialize state space model
%Input
u1 = 2 + sin(2*pi*t);
u2 = 2 + cos(2*pi*t);
u = [u1;u2];
%System
x = zeros(2,iter);
y = zeros(2,iter);
%Luenberger Observer
xo = zeros(2,iter);
yo = zeros(2,iter);
%Unknown Input Observer
z = zeros(2,iter);
xz = zeros(2,iter);
yz = zeros(2,iter);
%Observation Error
eo = zeros(2,iter);
ez = zeros(2,iter);

%White noise generation
w = randn(1,iter);
w = w/std(w);
w = w - mean(w);        
v = 0.5;
sigma = 0.1;
w = v + sigma * w;

%% Iterative Observation
for k = 1:iter    
    if k >= 150
        T(1,:) = 0;
        T(2,:) = 0.5*T(2,:);
    end
    
    %System
    y(:,k) = C * x(:,k);
    x(:,k+1) = A * x(:,k) + B * u(:,k) + E * w(k);
    
    %Luenberger Observer
    yo(:,k) = C * xo(:,k);
    xo(:,k+1) = (A - L * C) * xo(:,k) + B * u(:,k) + L * y(:,k);
    
    %Unknown Input Observer
    xz(:,k) = z(:,k) + H * y(:,k);
    z(:,k+1) = N * z(:,k) + T * u(:,k) + K * y(:,k);
    yz(:,k) = C * xz(:,k);
    
    %Residual Error
    eo(:,k) = y(:,k) - yo(:,k);
    ez(:,k) = y(:,k) - yz(:,k);
end

%% Plotting
figure(1)
plot(t,y,'--','linewidth',2);hold on
plot(t,yo,'-o',t,yz);hold on
line([1.5 1.5],[-1 9],'linestyle','--')
xlabel('t/s')
ylabel('Output')
legend('y_1','y_2','y_{1LBG}','y_{2LBG}','y_{1UIO}','y_{2UIO}')
figure(2)
plot(t,eo,'-o',t,ez);hold on
line([1.5 1.5],[-20 10],'linestyle','--')
xlabel('t/s')
ylabel('Residual Error')
legend('e_{1LBG}','e_{2LBG}','e_{1UIO}','e_{2UIO}')
figure(3)
plot(t,x(:,1:iter),'--','linewidth',2);hold on
plot(t,xo(:,1:iter),'-o',t,xz(:,1:iter));
xlabel('t/s')
ylabel('State Variables')
legend('x_1','x_2','x_{1LBG}','x_{2LBG}','x_{1UIO}','x_{2UIO}')
line([1.5 1.5],[-10 20],'linestyle','--')
figure(4)
subplot(211)
plot(t,ez(1,:))
line([1.5 1.5],[-20 10],'linestyle','--')
xlabel('t/s')
ylabel('Residual error')
legend('r_{1UIO}')
subplot(212)
plot(t,ez(2,:))
line([1.5 1.5],[-2 10],'linestyle','--')
xlabel('t/s')
ylabel('Residual error')
legend('r_{2UIO}')