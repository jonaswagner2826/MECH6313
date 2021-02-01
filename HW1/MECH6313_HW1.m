%% MECH6313 - HW 1
clear
close all

pblm1 = false;
pblm2 = true;

if pblm1
%% Problem 1
% System Def
delta = 0.05;
alpha = 0.4;
omega_t = 1.3;
sys1 = nlsys(@duff_eq);

% Simulation Setup
x_0 = [0;0];

N = 1e4;
t_step = 0.01;
t_max = N * t_step - t_step;
T = reshape(0:t_step:t_max,N,1);
U = alpha * cos(omega_t * T);
SYS1 = nlsim(sys1,U,T,x_0);

% Phase Plot
fig = SYS1.phasePlot(1,2,'Problem 1 - Phase Plot (Relaxed System)');
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm1_phase.png'))

% Phase Plot Comparrison
N = 5e3;
fig = figure('position',[0,0,1000,1000]);
for i = 1:4
    x_0 = randn(2,1);
    SYS(i) = nlsim(sys1,U,T,x_0);
    ax = subplot(2,2,i);
    SYS(i).phasePlot(1,2,x_0,fig,ax);
end
sgtitle('Problem 1 - Duff Equation Phase Portrait Comparrison')
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm1_phase_comparrision.png'))

% Vs Time Plot
fig = SYS1.plot(-1,0,0);
sgtitle('Problem 1 - Duff Equation Time Simulation')
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm1_vs_time.png'))

end

if pblm2
%% Problem 2
% System Def
sys2 = nlsys(@van_der_pol);

% Simulation Setup
x_0 = [0.8; -0.2];

N = 5e3;
t_step = 0.01;
t_max = N * t_step - t_step;
T = reshape(0:t_step:t_max,N,1);
U = 0 * T;
SYS2 = nlsim(sys2,U,T,x_0);

% Phase Plot
fig = SYS2.phasePlot(1,2,'Problem 2 - Phase Plot');
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm2_phase.png'))

% Phase Plot Comparrison
fig = figure('position',[0,0,1000,1000]);
X_0 = [ -0.5,0.8, -1.5, 3;
        0.5, -0.5, 2.7, -1.9];
for i = 1:4
    SYS(i) = nlsim(sys2,U,T,X_0(:,i));
    ax = subplot(2,2,i);
    SYS(i).phasePlot(1,2,X_0(:,i),fig,ax);
end
sgtitle('Problem 2 - Van Der Pol Phase Portraits')
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm2_phase_comparrision.png'))

% Vs Time Plot
x_0 = X_0(1);
fig = SYS2.plot(-1,0,0);
sgtitle('Problem 2 - Van Der Pol Time Simulation')
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm2_vs_time.png'))


% Negative Equivelent
sys2_neg= nlsys(@van_der_pol_neg);

% Negative Stability
x_0 = [0;0];
u_0 = 0;
sys2_neg_lin = sys2_neg.ss(x_0,u_0)
eig_A = eig(sys2_neg_lin.A)
is_stable = isstable(sys2_neg_lin)



%Negative Sim
x_0 = [0.8; -0.2];
SYS2_neg = nlsim(sys2_neg,U,T,x_0)


% Phase Plot Comparrison
fig = figure('position',[0,0,1000,1000]);
%X_0(i) is from last set of plot
for i = 1:4
    SYS_neg(i) = nlsim(sys2_neg,U,T,X_0(:,i));
    ax = subplot(2,2,i);
    SYS_neg(i).phasePlot(1,2,X_0(:,i),fig,ax);
end
sgtitle('Problem 2 - Negative Van Der Pol Phase Portraits')
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm2_phase_comparrision_neg.png'))

% Vs Time Plot
x_0 = X_0(1);
fig = SYS2_neg.plot(-1,0,0);
sgtitle('Problem 2 - Negative Van Der Pol Time Simulation')
saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm2_vs_time_neg.png'))


end


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

function y = van_der_pol(x,u,a)
    % VAD_DER_POL nonlin function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        a = 1;
    end
    
    % Array sizes
    n = 2; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = x(2);
    y(2,1) = - a * (x(1)^2 -1) * x(2) - x(1) + u;
    
    
    if nargin ==0
        y = [n;p];
    end
end

function y = van_der_pol_neg(x,u,a)
    % VAD_DER_POL nonlin function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        a = -1;
    end
    
    % Array sizes
    n = 2; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = x(2);
    y(2,1) = - a * (x(1)^2 -1) * x(2) - x(1) + u;
    
    
    if nargin ==0
        y = [n;p];
    end
end