% MECH 6313 - Exam

clear
close all

pblm1 = true;
pblm2 = true;
pblm3 = false;
pblm4 = false;
pblm5 = false;
pblm6 = false;

if pblm1
%% Problem 1

satVal = 2;

tau = 0.25;
mu = 0.1;
% Phase Plot 1
figure('position',[0,0,1200,1000])
xmax = 5;
ymax = 5;
xmin = -5;
ymin = -5;
xstep = 0.25;
ystep = xstep;

[X,Y] = meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
DX = max(min((X - (X.^3)/3 - Y)/tau, satVal), -satVal);
DY = max(min(X + mu, satVal), -satVal);

quiver(X,Y,DX,DY)
title('Phase Portrait \mu < \mu_c')
hold on
% x = [xmin:xstep:xmax];
% y = max(min(1/2*tan(pi/2 * x), xmax), xmin);
% plot(x,y, 'LineWidth', 2)
% plot(y,x, 'LineWidth', 2)
saveas(gcf,[pwd,'\Exam\fig\pblm1_phaseplot_mu01.png'])



% tau = 0.1;
mu = 2;
% Phase Plot 2
figure('position',[0,0,1200,1000])
% xmax = 1;
% ymax = xmax;
% xmin = -1.5*xmax;
% ymin = xmin;
% xstep = 0.5;
% ystep = xstep;

[X,Y] = meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
DX = max(min((X - (X.^3)/3 - Y)/tau, satVal), -satVal);
DY = max(min(X + mu, satVal), -satVal);

quiver(X,Y,DX,DY)
title('Phase Portrait \mu > \mu_c')
hold on
% x = [xmin:xstep:xmax];
% y = max(min(1/2*tan(pi/2 * x), xmax), xmin);
% plot(x,y, 'LineWidth', 2)
% plot(y,x, 'LineWidth', 2)
saveas(gcf,[pwd,'\Exam\fig\pblm1_phaseplot_mu2.png'])


end

if pblm2
%% Problem 2
solveEqPnt = false;
phasePlt = true;
linSysCalc = false;


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
figure('position',[0,0,1200,1000])
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
saveas(gcf,[pwd,'\Exam\fig\pblm2_phaseplot.png'])


% Phase Plot 2
figure('position',[0,0,1200,1000])
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
hold on
% Region of Attraction
viscircles([[0.5,0.5;-0.5,-0.5]],[0.5,0.5])
x=[-1,-1,1,1,-1,1];
y=[1,-1,-1,1,1,-1];
plot(x,y,'c','LineWidth',4)
% Eq-Points
scatter([-0.5,0.5],[-0.5,0.5],100,'rp')
scatter([0],[0],100,'kp')
saveas(gcf,[pwd,'\Exam\fig\pblm2_phaseplot_origin.png'])


% Phase Plot 3
figure('position',[0,0,1200,1000])
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
saveas(gcf,[pwd,'\Exam\fig\pblm2_phaseplot_zoomOut.png'])
end



if linSysCalc
% -------------------------------------------------
% Linearized System Calc
sys2a = nlsys(@pblm2a)
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

syms x1 x2
eq1 = 0 == -x1 + x2^3 + 1;
eq2 = 0 == -4*x1^2 + 3*x2;
solve([eq1,eq2],[x1,x2])

sys5 = nlsys(@pblm5a)

% Phase Plot 2
figure()
xmax = 5;
ymax = xmax;
xmin = -xmax;
ymin = -ymax;
xstep = 0.1;
ystep = xstep;

[X,Y] = meshgrid(xmin:xstep:xmax,ymin:ystep:ymax);
DX = -X + Y^3 + 1%max(min(, 1), -1);
DY = -4*X^2 + 3*Y%max(min(, 1), -1);

quiver(X,Y,DX,DY)


end

if pblm6
%% Problem 6
sys6 = nlsys(@pblm6a)

linsys6 = sys6.linearize([0;0])


end


%% Local Functions
function y = pblm1a(x,u)
    % pblm1c function
    arguments
        x (2,:) = [0; 0];
        u (1,:) = 0;
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    
    % Parameters
    tau = 0.1;
    mu = 0.9;
    
    
    % State Upadate Eqs
    y(1,1) = (x(1) - (x(1)^3)/3 - x(2))/tau;
    y(2,1) = x(1) + mu;
    
    if nargin == 0
        y = [n;p];
    end 
end

function y = pblm1b(x,u)
    % pblm1c function
    arguments
        x (2,:) = [0; 0];
        u (1,:) = 0;
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    % Parameters
    tau = 0.1;
    mu = 1.1;
    
    
    % State Upadate Eqs
    y(1,1) = (x(1) - (x(1)^3)/3 - x(2))/tau;
    y(2,1) = x(1) + mu;
    
    if nargin == 0
        y = [n;p];
    end 
end


function y = pblm2b(x,u)
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

function y = pblm5a(x,u)
    % pblm1c function
    arguments
        x (2,:) = [0; 0];
        u (1,:) = 0;
    end
    
    % Array Sizes
    n = 2;
    p = 1;
    
    
    % State Upadate Eqs
    y(1,1) = -x(1) + x(2)^3 + 1;
    y(2,1) = -4*x(1)^2 + 3*x(2);
    
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