
%% creates a system of ODEs
function dy = example(t,y,a)

dy = zeros(2,1);  % state: a column vector

% Hopf bifurcation example - 2nd order system
dy(1) = y(1)*(a - y(1)^2 - y(2)^2) - y(2);
dy(2) = y(2)*(a - y(1)^2 - y(2)^2) + y(1);

