%% MECH6313 - HW 3
clear
close all

pblm1 = true;
pblm2 = false;
pblm3 = false;




if pblm1
%% Problem 1
syms x_1 x_2
x_1_dot = -x_1 + x_2
x_2_dot = (x_1^2)/(1 + x_1^2) - 0.5 * x_2

% part a
syms x_1_bar x_2_bar

x_1_bar_dot = subs(x_1_dot, [x_1, x_2], [x_1_bar + 1, x_2_bar + 1]);
x_2_bar_dot = subs(x_2_dot, [x_1, x_2], [x_1_bar + 1, x_2_bar + 1]);

x_bar = [x_1_bar; x_2_bar];
x_bar_dot = [x_1_bar_dot; x_2_bar_dot]

% part b
% Linearize
A_sym = jacobian(x_bar_dot, x_bar)
A = subs(A_sym, [x_1_bar, x_2_bar], [0,0])

[T1, eig_A] = eig(A)

% Transform
syms y_sym z_sym
x_bar_sub = [y_sym, z_sym] * T1;
y_dot = subs(x_1_bar_dot, [x_1_bar, x_2_bar], x_bar_sub);
z_dot = subs(x_2_bar_dot, [x_1_bar, x_2_bar], x_bar_sub);

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
w_dot = A2 * h_sym + subs(g2,z_sym,h_sym) - dh_sym * (A1 * y_sym + subs(g1,z_sym,h_sym));

% Taylor's Series approximation of manifold
syms h2 h3
h = h2 * y_sym^2;% + h3 * y_sym^3;
dh = diff(h,y_sym);

% Attempting to solve for the h2 value...
w_dot = expand(subs(w_dot, [h_sym, dh_sym], [h, dh]));
'w_dot'
pretty(w_dot)

w_dot_soln = expand(w_dot);
w_dot_soln = subs(w_dot_soln, y_sym^4, 0);
w_dot_soln = subs(w_dot_soln, y_sym^3, 0);
% w_dot_soln = subs(w_dot_soln, y_sym^2, 0);
% w_dot_soln = subs(w_dot_soln, y_sym, 0);
'w_dot_soln'
pretty(w_dot_soln)


end

%% Local Functions
% function dx = pblm1a(t, x, parms)
%     % pblm1a function
%     arguments
%         t (1,1) = 0;
%         x (2,1) = [0; 0];
%         parms = false;
%     end
%     
%     if parms == false
%         alpha = 1;
%     else
%         alpha = parms(1);
%     end
% 
%     % State Upadate Eqs
%     dx(1,1) = alpha * x(1) + x(2);
%     dx(2,1) = - x(1) + alpha*x(2) - x(1)^2 * x(2);
% end
% 
% function y = pblm1b(t,x,parms)
%     % pblm1b function
%     arguments
%         t (1,1) = 0;
%         x (2,1) = [0; 0];
%         parms = false;
%     end
%     
%     if parms == false
%         alpha = 1;
%     else
%         alpha = parms(1);
%     end
%     
%     % State Upadate Eqs
%     y(1,1) = alpha * x(1) + x(2) - x(1)^3;
%     y(2,1) = - x(1) + alpha*x(2) + 2 *x(2)^3;
% end
% 
% function y = pblm1c(t,x,parms)
%     % pblm1c function
%     arguments
%         t (1,1) = 0;
%         x (2,1) = [0; 0];
%         parms = false;
%     end
%     
%     if parms == false
%         alpha = 1;
%     else
%         alpha = parms(1);
%     end
%     % State Upadate Eqs
%     y(1,1) = alpha * x(1) + x(2) - x(1)^2;
%     y(2,1) = - x(1) + alpha*x(2) + 2 * x(1)^2;
% end
% 
% function y = pblm2a(x,u)
%     % pblm2 function
%     arguments
%         x (2,1) = [0; 0];
%         u (1,1) = 0;
%     end
%     
%     % Array Sizes
%     n = 2;
%     p = 1;
%     
%     % State Upadate Eqs
%     y(1,1) = x(2) + x(1) * x(2)^2;
%     y(2,1) = - x(1) + x(1)^2 * x(2);
%     
%     if nargin == 0
%         y = [n;p];
%     end
% end