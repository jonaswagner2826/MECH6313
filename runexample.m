
%% Explains how to run an ode system
% simulation parameters
a = 1;

% initial condition
x0 = [10; 15]; 

% simulation interval
Tspan = 0:0.01:10;

options = odeset('RelTol',1e-6,'AbsTol',1e-9);

tic
[T,Y] = ode45(@example,Tspan,x0,options,a);
toc

% plot states in the phase plane
figure()
plot(Y(:,1),Y(:,2),'*')

% plot 1st state vs time
figure()
plot(T,Y(:,1))

% plot 2nd state vs time
figure()
plot(T,Y(:,2))


