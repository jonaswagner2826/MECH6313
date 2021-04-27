% MECH 6313 - HW6
% Jonas Wagner
% 2021-04-27
% 

clear
close all

pblm1 = false;
pblm2 = true;


if pblm1
%% Problem 1

G1 = tf([5 3 1], [1 2 1])
isPassive(G1)

G2 = tf([1 1 5 0.1], [1 2 3 4])
isPassive(G2)

Gp = G1 + G2
isPassive(Gp)

Gs = G1 * G2
isPassive(Gs)

end


if pblm2
%% Problem 2
pblm2a = false;
pblm2b = true;
pblm2c = true;

if pblm2a
% Part a
syms omega a b lambda
assume(a,'real'); assume(a > 0)
assume(b,'real'); assume(b > 0)
assume(omega,'real')
assume(lambda, 'real')

num = j*omega + lambda;
den = omega^2 + j*a*omega + b;

H_sym = num/den;
disp('H(s) = ')
pretty(H_sym)

H_real = real(H_sym);
disp('H_real = ')
pretty(H_real)

H_imag = imag(H_sym);
disp('H_imag = ')
pretty(imag(H_sym))
end

if pblm2b
% Part b
lambda1 = 1;
lambda2 = -1;

H1 = tf([1 lambda1], [1 1 1])
isPassive(H1)
figure
nyquist(H1)
title('Nyquist Plot for H_1(s)')
saveas(gcf, [pwd, '\Homework\HW6\fig\pblm2_H1.png'])

H2 = tf([1 lambda2], [1 1 1])
isPassive(H2)
figure
nyquist(H2)
title('Nyquist Plot for H_2(s)')
saveas(gcf, [pwd, '\Homework\HW6\fig\pblm2_H2.png'])
end

if pblm2c
H2_sys = ss([0, 1; -1, -1], [0; 1], [-1 1], 0)
tf (H2_sys)
end


end
