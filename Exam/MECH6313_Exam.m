% MECH 6313 - Exam

clear
close all

pblm1 = false;
pblm2 = false;
pblm3 = false;
pblm4 = false;
pblm5 = false;
pblm6 = true;

if pblm1
%% Problem 1
end

if pblm2
%% Problem 2
solveEqPnt = true;
phasePlt = true;
linSysCalc = true;


if solveEqPnt
% --------------------------------------------
% Equalibrium Points
syms x1 x2
eq1 = 0 == -1/2 * tan(pi*x1/2) + x2;
eq2 = 0 == x1 -1/2 * tan(pi*x1/2);


[x1_eq, x2_eq] = vpasolve([eq1, eq2], [x1,x2]);

eq3 = x1 == 1/2 * tan(pi*x1/2);

x3_eq = solve(eq3, x1);
end

if phasePlt
% ----------------------------------------------
% Phase Plot 1
figure()
xmax = 2.5;
ymax = xmax;
xmin = -xmax;
ymin = -ymax;
xstep = 0.1;
ystep = xstep;

[X,Y] = meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
DX = max(min(-1/2 * tan(pi*X/2) + Y, 1), -1);
DY = max(min(-1/2 * tan(pi*Y/2) + X, 1), -1);

quiver(X,Y,DX,DY)
title('Phase Portrait')
hold on
x = [xmin:xstep:xmax];
y = max(min(1/2*tan(pi/2 * x), xmax), xmin);
plot(x,y, 'LineWidth', 2)
plot(y,x, 'LineWidth', 2)


% Phase Plot 2
figure()
xmax = 1;
ymax = xmax;
xmin = -xmax;
ymin = -ymax;
xstep = 0.1;
ystep = xstep;

[X,Y] = meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
DX = max(min(-1/2 * tan(pi*X/2) + Y, 1), -1);
DY = max(min(-1/2 * tan(pi*Y/2) + X, 1), -1);

quiver(X,Y,DX,DY)
title('Phase Portrait (Zoomed-In)')
% hold on
% x = [xmin:xstep:xmax];
% y = max(min(1/2*tan(pi/2 * x), xmax), xmin);
% plot(x,y, 'LineWidth', 2)
% plot(y,x, 'LineWidth', 2)


% Phase Plot 3
figure()
xmax = 4;
ymax = xmax;
xmin = -xmax;
ymin = -ymax;
xstep = 0.1;
ystep = xstep;

[X,Y] = meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
DX = max(min(-1/2 * tan(pi*X/2) + Y, 1), -1);
DY = max(min(-1/2 * tan(pi*Y/2) + X, 1), -1);

quiver(X,Y,DX,DY)
title('Phase Portrait (Zoomed-Out)')
hold on
x = [xmin:xstep:xmax];
y = max(min(1/2*tan(pi/2 * x), xmax), xmin);
plot(x,y, 'LineWidth', 2)
plot(y,x, 'LineWidth', 2)



% U = X;
% V = 0.5*Y;

% 
% 
% sys def
sys2a = nlsys(@pblm2a)
% 
% % Simulation Setup
% x_0 = [-0.1;-0.05];
% 
% N = 500;
% t_step = 0.01;
% t_max = N * t_step - t_step;
% T = reshape(0:t_step:t_max,N,1);
% U = 0*T;
% SYS2 = nlsim(sys2a,U,T,x_0);

% Phase Plot
% fig = SYS2.phasePlot(1,2,'Problem 1 - Phase Plot (Relaxed System)');
end

if linSysCalc
% -------------------------------------------------
% Linearized System Calc
syms x1 x2
linsys2a_sym = sys2a.linearize([x1;x2])
linsys2_0 = sys2a.linearize([0;0])
eig(linsys2_0)
linsys2_p05 = sys2a.linearize([0.5;0.5])
eig(linsys2_p05)
linsys2_n05 = sys2a.linearize([-0.5;-0.5])
eig(linsys2_n05)
linsys2_p1p1 = sys2a.linearize([1;1])
linsys2_n1n1 = sys2a.linearize([-1;-1])
linsys2_p1n1 = sys2a.linearize([1;-1])
linsys2_p125n125 = sys2a.linearize([1.25;-1.25])
eig(linsys2_p125n125)
linsys2_n1p1 = sys2a.linearize([-1;1])
linsys2_n125p125 = sys2a.linearize([-1.25;1.25])
eig(linsys2_n125p125)
end



end

if pblm3
%% Problem 3
end

if pblm4
%% Problem 4
end

if pblm5
%% Problem 5
end

if pblm6
%% Problem 6
sys6 = nlsys(@pblm6a)

linsys6 = sys6.linearize([0;0])


end


%% Local Functions
function y = pblm2a(x,u)
    % pblm1c function
    arguments
        x (2,:) = [0; 0];
        u (1,:) = 0;
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    
    % State Upadate Eqs
    y(1,1) = -1/2 * tan(pi*x(1)/2) + x(2);
    y(2,1) = x(1) -1/2 * tan(pi*x(2)/2);
    
    if nargin == 0
        y = [n;p];
    end 
end

function y = pblm6a(x,u)
    % pblm1c function
    arguments
        x (2,:) = [0; 0];
        u (1,:) = 0;
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    
    % State Upadate Eqs
    y(1,1) = x(2);
    y(2,1) = -(sin(x(1)) + 2) * (x(1) + x(2));
    
    if nargin == 0
        y = [n;p];
    end 
end