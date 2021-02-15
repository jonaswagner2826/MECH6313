%% MECH6313 - HW 2
clear
close all

pblm1 = false;
pblm2 = true;


if pblm1
%% Problem 1
parta = false;
partb = false;
partc = true;

if parta
%% Problem 1a

% System Def
% alpha < 0
alpha = -1;
sys1a_neg = nlsys(@pblm1a,'empty',0,-1,0,0,alpha);
% alpha = 0
alpha = 0;
sys1a_zero = nlsys(@pblm1a,'empty',0,-1,0,0,alpha);
% alpha > 0
alpha = 1;
sys1a_pos = nlsys(@pblm1a,'empty',0,-1,0,0,alpha);
%Sys Array
sysArray = [sys1a_neg, sys1a_zero, sys1a_pos];

% Simulation Setup
N = 5e3;
t_step = 0.01;
t_max = N * t_step - t_step;
T = reshape(0:t_step:t_max,N,1);
U = 0 * T;

% Phase Plot Comparrison
fig = figure('position',[0,0,1500,500]);
X_0 = [ -0.5,0.8, -1.5, 3;
        0.5, -0.5, 2.7, -1.9];

% Sim Phase Plots
simNum = 1;
for i = 1:3
    sys = sysArray(i);
    ax = subplot(1,3,i);
    for j = 1:4
        SYS(simNum) = nlsim(sys,U,T,X_0(:,j));
        SYS(simNum).phasePlot(1,2,'',fig,ax);
        hold on
        simNum = simNum + 1;
    end
end

subplot(1,3,1)
title('$\alpha$ = -1','interpreter','latex')
subplot(1,3,2)
title('$\alpha$ = 0','interpreter','latex')
subplot(1,3,3)
title('$\alpha$ = 1','interpreter','latex')

sgtitle('Problem 1a - Phase Portrait For Varyuing Alpha')
saveas(fig,fullfile([pwd '\\' 'HW2' '\\' 'fig'],'pblm1a.png'))

end

if partb
%% Problem 1b

% System Def
% alpha < 0
alpha = -1;
sys1b_neg = nlsys(@pblm1b,'empty',0,-1,0,0,alpha);
% alpha = 0
alpha = 0;
sys1b_zero = nlsys(@pblm1b,'empty',0,-1,0,0,alpha);
% alpha > 0
alpha = 1;
sys1b_pos = nlsys(@pblm1b,'empty',0,-1,0,0,alpha);
%Sys Array
sysArray = [sys1b_neg, sys1b_zero, sys1b_pos];

% Simulation Setup
N = 5e3;
t_step = 0.01;
t_max = N * t_step - t_step;
T = reshape(0:t_step:t_max,N,1);
U = 0 * T;

% Phase Plot Comparrison
fig = figure('position',[0,0,1500,500]);
X_0 = [ -0.5,0.8, -0.5, 0.3;
        -0.2, -0.5, 0.2, -0.9];

% Sim Phase Plots
simNum = 1;
for i = 1:3
    sys = sysArray(i);
    ax = subplot(1,3,i);
    for j = 1:4
        SYS(simNum) = nlsim(sys,U,T,X_0(:,j));
        SYS(simNum).phasePlot(1,2,'',fig,ax);
        hold on
        simNum = simNum + 1;
    end
end

subplot(1,3,1)
title('$\alpha$ = -1','interpreter','latex')
subplot(1,3,2)
title('$\alpha$ = 0','interpreter','latex')
subplot(1,3,3)
title('$\alpha$ = 1','interpreter','latex')

sgtitle('Problem 1b - Phase Portrait For Varyuing Alpha')
saveas(fig,fullfile([pwd '\\' 'HW2' '\\' 'fig'],'pblm1b.png'))

end


if partc
%% Problem 1c

% System Def
% alpha < 0
alpha = -1;
sys1c_neg = nlsys(@pblm1c,'empty',0,-1,0,0,alpha);
% alpha = 0
alpha = 0;
sys1c_zero = nlsys(@pblm1c,'empty',0,-1,0,0,alpha);
% alpha > 0
alpha = 1;
sys1c_pos = nlsys(@pblm1c,'empty',0,-1,0,0,alpha);
%Sys Array
sysArray = [sys1c_neg, sys1c_zero, sys1c_pos];

% Simulation Setup
N = 5e3;
t_step = 0.01;
t_max = N * t_step - t_step;
T = reshape(0:t_step:t_max,N,1);
U = 0 * T;

% Phase Plot Comparrison
fig = figure('position',[0,0,1500,500]);
X_0 = [ -0.5,0.8, -0.5, 0.3;
        -0.2, -0.5, 0.2, -0.9];

% Sim Phase Plots
simNum = 1;
for i = 1:3
    sys = sysArray(i);
    ax = subplot(1,3,i);
    for j = 1:4
        SYS(simNum) = nlsim(sys,U,T,X_0(:,j));
        SYS(simNum).phasePlot(1,2,'',fig,ax);
        hold on
        simNum = simNum + 1;
    end
end

subplot(1,3,1)
title('$\alpha$ = -1','interpreter','latex')
subplot(1,3,2)
title('$\alpha$ = 0','interpreter','latex')
subplot(1,3,3)
title('$\alpha$ = 1','interpreter','latex')

sgtitle('Problem 1c - Phase Portrait For Varyuing Alpha')
saveas(fig,fullfile([pwd '\\' 'HW2' '\\' 'fig'],'pblm1c.png'))
end
    
end

if pblm2
%% Problem 2
parta = true;

if parta
%% Problem 2a
% sys def
sys2a = nlsys(@pblm2a);

syms x1 x2
linsys2a_sym = sys2a.linearize([x1;x2])
linsys2a = sys2a.linearize([0;0])

end

end





%% Local Functions
function y = pblm1a(x,u,parms)
    % pblm1a function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        alpha = 1;
    else
        alpha = parms(1);
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    % State Upadate Eqs
    y(1,1) = alpha * x(1) + x(2);
    y(2,1) = - x(1) + alpha*x(2) - x(1)^2 * x(2);
    
    if nargin == 0
        y = [n;p];
    end
end

function y = pblm1b(x,u,parms)
    % pblm1a function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        alpha = 1;
    else
        alpha = parms(1);
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    % State Upadate Eqs
    y(1,1) = alpha * x(1) + x(2) - x(1)^3;
    y(2,1) = - x(1) + alpha*x(2) + 2 *x(2)^2;
    
    if nargin == 0
        y = [n;p];
    end
end

function y = pblm1c(x,u,parms)
    % pblm1a function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
        parms = false
    end
    
    if parms == false
        alpha = 1;
    else
        alpha = parms(1);
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    % State Upadate Eqs
    y(1,1) = alpha * x(1) + x(2) - x(1)^2;
    y(2,1) = - x(1) + alpha*x(2) - x(1)^2 * x(1)^2;
    
    if nargin == 0
        y = [n;p];
    end
end

function y = pblm2a(x,u)
    % pblm2 function
    arguments
        x (2,1) = [0; 0];
        u (1,1) = 0;
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    % State Upadate Eqs
    y(1,1) = x(2) + x(1) * x(2)^2;
    y(2,1) = - x(1) + x(1)^2 * x(2);
    
    if nargin == 0
        y = [n;p];
    end
end
