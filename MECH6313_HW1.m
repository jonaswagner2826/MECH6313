%% MECH6313 - HW 1
clear
close all



%% Problem 1
% System Def
delta = 0.05;
alpha = 0.4;
omega_t = 1.3;
sys1 = nlsys(@duff_eq,delta);

% Simulation Setup
x_0 = [0;0];

N = 5000;
t_step = 0.01;
t_max = N * t_step - t_step;
T = reshape(0:t_step:t_max,N,1);
U = alpha * cos(omega_t * T);
SYS1 = nlsim(sys1,U,T,x_0);

% Phase Plot
SYS1.phasePlot(1,2,'Problem 1 - Phase Plot (Relaxed System)');

% Phase Plot Comparrison
n = 5;
fig = figure('position',[200,300,2000,500]);
for i = 1:n
    x_0 = randn(2,1);
    SYS(i) = nlsim(sys1,U,T,x_0);
    ax = subplot(1,n,i);
    SYS(i).phasePlot(1,2,x_0,fig,ax);
end

% Vs Time Plot
SYS1.plot()


%% Local Functions
function y = duff_eq(x,u,delta)
    % DUFF_EQ nonlin function with caotic behavior
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        delta = 0.05;
    end
    
    % Array sizes
    n = 2; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = x(2);
    y(2,1) = - delta * x(2) + x(1) - x(1)^3 + u;
    
    
    if nargin ==0
        y = [n;p];
    end
end