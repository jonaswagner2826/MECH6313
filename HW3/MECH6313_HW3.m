%% MECH6313 - HW 3
clear
close all

pblm1 = true;
pblm2 = false;
pblm3 = false;
pblm4 = false;
pblm5 = false;


if pblm1
%% Problem 1
syms x_1 x_2
x_1_dot = -x_1 + x_2
x_2_dot = (x_1^2)/(1 + x_1^2) - 0.5 * x_2

% part a
syms x_1_bar x_2_bar

x_1_bar_dot = subs(x_1_dot, [x_1, x_2], [x_1_bar + 1, x_2_bar + 1])
x_2_bar_dot = subs(x_2_dot, [x_1, x_2], [x_1_bar + 1, x_2_bar + 1])

x_bar = [x_1_bar; x_2_bar]
x_bar_dot = [x_1_bar_dot; x_2_bar_dot]

% part b
% Linearize
A_sym = jacobian(x_bar_dot, x_bar)
A = subs(A_sym, [x_1_bar, x_2_bar], [0,0])

[T1, eig_A] = eig(A)

% Transform
syms y_sym z_sym
assume(y_sym,'real')
assume(z_sym,'real')
x_bar_sub = T1 * [y_sym; z_sym];
y_dot = subs(x_1_bar_dot, [x_1_bar, x_2_bar], x_bar_sub');
z_dot = subs(x_2_bar_dot, [x_1_bar, x_2_bar], x_bar_sub');

% G function Definitions
g1 = y_dot;
g2 = z_dot + 3/2 * z_sym;
'g1'
pretty(g1)
'g2'
pretty(g2)

% Coefficents of eigenvalue matrix
A1 = eig_A(1,1);
A2 = eig_A(2,2);

% w_dot substitution (from definition equation)
syms h_sym dh_sym
w_dot = A2 * h_sym + subs(g2,z_sym,h_sym) == dh_sym * (A1 * y_sym + subs(g1,z_sym,h_sym));

% Taylor's Series approximation of manifold
syms h2
h = h2 * y_sym^2;
dh = diff(h,y_sym);

% Attempting to solve for the h2 value...
w_dot = (subs(w_dot, [h_sym, dh_sym], [h, dh]));
'w_dot'
pretty(w_dot)

solve(w_dot,h2)

% w_dot_soln = expand(w_dot);
% w_dot_soln = subs(w_dot_soln, y_sym^4, 0);
% w_dot_soln = subs(w_dot_soln, y_sym^3, 0);
% syms y2
% w_dot_soln = subs(w_dot_soln, y_sym^2, y2);
% w_dot_soln = subs(w_dot_soln, y_sym, 0);
% 'w_dot_soln'
% pretty(w_dot_soln)
% 
% syms h2y2
% w_dot_soln = subs(w_dot_soln, h2*y2, h2y2)
% solve(w_dot_soln == 0,h2y2)


end



if pblm4
%% Problem 4
syms x1 x2
x1_dot = -x1 + (2*x2)/(1 + x2^2);
x2_dot = -x2 + (2*x1)/(1 + x1^2);
f = [x1_dot; x2_dot];
'f'
pretty(f)
df = jacobian(f);
'jacobian'
pretty(df)




end


if pblm5
%% Problem 5
syms x1 x2 a b c
x1_dot = atan(a * x1) - x1 * x2;
x2_dot = b * x1^2 - c * x2;
x_dot = [x1_dot; x2_dot];
'x_dot'
pretty(x_dot)

x = [x1; x2];
'x'
pretty(x)

mu = [a; b; c];
mu_bar = [1; 0; 1];
'mu'
pretty(mu)


A_tau = jacobian(x_dot, x)
B_tau = jacobian(x_dot, mu)


sys_func = @pblm5_func;

T = [0,10];
x_0 = [1,-1, 0,0,0,0,0,0]';

[t,y] = ode45(@(t,y) sys_func(t,y,mu_bar,A_tau,B_tau),T,x_0);

y_states = y(:,[1,2]);
y_a = y(:,[3,6]);
y_b = y(:,[4,7]);
y_c = y(:,[5,8]);

fig = figure('position',[0,0,1500,1200]);
subplot(3,3,[1:3])
plot(t,y_states)
title('States')

titles = ["f1 - a", "f1 - b", "f1 - c", "f2 - a", "f2 - b", "f2 - c"];
for i = 3:8
    subplot(3,3,i+1)
    plot(t,y(:,i))
    title(titles(i-2))
end

saveas(fig,fullfile([pwd '\\' 'HW3' '\\' 'fig'],'pblm5.png'))

end





%% Local Functions
function dx = pblm5_func(t, x, parms, A, B)
    % pblm5 function
    arguments
        t (1,1) = 0;
        x (8,1) = zeros(8,1); %state and 6 sensitivities
        parms = false;
        A = 0;
        B = 0;
    end
    
    if parms == false
        a = 1;
        b = 0;
        c = -1;
    else
        a = parms(1);
        b = parms(2);
        c = parms(3);
    end
    
    % Variable Decode
    x1 = x(1);
    x2 = x(2);
    S = zeros(2,3);%[x(3), x(4), x(5); x(6), x(7), x(8)];

    % State Upadate Eqs
    x1_dot = atan(a * x1) - x1 * x2;
    x2_dot = b * x1^2 - c * x2;
    S_dot = subs(A * S + B);
    
    % Variable Encode
    dx = x;
    dx(1) = x1_dot;
    dx(2) = x2_dot;
    dx(3) = S_dot(1,1);
    dx(4) = S_dot(1,2);
    dx(5) = S_dot(1,3);
    dx(6) = S_dot(2,1);
    dx(7) = S_dot(2,2);
    dx(8) = S_dot(2,3);    
end

