%% MECH6313 - HW 1
clear
close all

pblm1 = false;
pblm2 = false;
pblm4 = true;

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
sys2_neg= nlsys(@van_der_pol,'empty',0,-1,0,0,-1);

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


if pblm4
%% Problem 4

syms f(x,r) g(x,r)

% Problem 3.4.2
sys342 = nlsys(@pblm342);


r = linspace(-2,5,20);
x = linspace(-4,4,20);
[r_c,fig] = sys342.bifurcationPlot(r,x);


xlimit = xlim;
plot([xlimit(1),r_c],[0,0],'g');
plot([r_c,xlimit(2)],[0,0],'g--');

x = linspace(r_c,xlimit(2));
y = 2 * sqrt(x-r_c);
plot(x,y,'g-')
plot(x,-y,'g-')
hold off

title('Problem 3.4.2')

saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm4_342.png'))


% Problem 3.4.4
sys344 = nlsys(@pblm344);
r = linspace(-5,3,20);
x = linspace(-4,4,20);
[r_c,fig] = sys344.bifurcationPlot(r,x);


xlimit = xlim;
plot([xlimit(1),r_c],[0,0],'g')
plot([r_c,xlimit(2)],[0,0],'g--')

x = linspace(xlimit(1),r_c);
y = 1 * sqrt(abs(x-r_c));
plot(x,y,'g--')
plot(x,-y,'g--')
hold off

title('Problem 3.4.4')

saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm4_344.png'))

% Problem 3.4.7
sys347 = nlsys(@pblm347);
r = linspace(-5,20,25);
x = linspace(-3,3,10);
[r_c,fig] = sys347.bifurcationPlot(r,x);


xlimit = xlim;
x = linspace(r_c,xlimit(2));
y = 0.35 * sqrt(x-r_c);
plot(x,y,'g--')
plot(x,-y,'g-')
hold off

title('Problem 3.4.7')

saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm4_347.png'))

% Problem 4.4.9
sys349 = nlsys(@pblm349);
r = linspace(-8,2,20);
x = linspace(-4,4,20);
[r_c,fig] = sys349.bifurcationPlot(r,x);


xlimit = xlim;
x = linspace(xlimit(1),-2);
y = 0 * x;
plot(x,y-1,'g--');
plot(x,y+1,'g--');
plot(x,y,'g-');

x = linspace(-2,r_c);
y = 0 * x;
plot(x,y-x-1,'g--');
plot(x,y+x+1,'g--');
plot(x,y,'g-');

x = linspace(r_c,xlimit(2));
y = 0*x;
plot(x,y,'g--');

hold off

title('Problem 3.4.9')

saveas(fig,fullfile([pwd '\\' 'HW1' '\\' 'fig'],'pblm4_349.png'))

end







%% Local Functions
function y = duff_eq(x,u,parms)
    % DUFF_EQ nonlin function with caotic behavior
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        delta = 0.05;
    else
        delta = parms(1);
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

function y = van_der_pol(x,u,parms)
    % VAD_DER_POL nonlin function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        a = 1;
    else
        a = parms(1);
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


function y = pblm342(x,u,parms)
    % VAD_DER_POL nonlin function
    arguments
        x (1,1) = 0;
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        r = 1;
    else
        r = parms(1);
    end
    
    % Array sizes
    n = 1; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = r * x(1) - sinh(x(1));
    
    
    if nargin ==0
        y = [n;p];
    end
end


function y = pblm344(x,u,parms)
    % VAD_DER_POL nonlin function
    arguments
        x (1,1) = 0;
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        r = 1;
    else
        r = parms(1);
    end
    
    % Array sizes
    n = 1; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = x(1) + (r*x(1))/(1+x(1)^2);
    
    
    if nargin ==0
        y = [n;p];
    end
end


function y = pblm347(x,u,parms)
    % VAD_DER_POL nonlin function
    arguments
        x (1,1) = 0;
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        r = 1;
    else
        r = parms(1);
    end
    
    % Array sizes
    n = 1; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = 5 - r * exp(-x(1)^2);
    
    
    if nargin ==0
        y = [n;p];
    end
end

function y = pblm349(x,u,parms)
    % VAD_DER_POL nonlin function
    arguments
        x (1,1) = 0;
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        r = 1;
    else
        r = parms(1);
    end
    
    % Array sizes
    n = 1; % Number of states
    p = 1; % Number of inputs

    % State Update Equations
    y(1,1) = x + tanh(r*x(1));
    
    
    if nargin ==0
        y = [n;p];
    end
end