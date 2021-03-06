% MECH 6313 - HW5
% Jonas Wagner
% 2021-04-08


% mu = 15; %Change this... or make it a global variable for looping...

%% Settings
% generateModel = false;
openModel = true;
simulateModel = true;
plotResults = true;

% Name of the simulink model
[cfolder,~,~] = fileparts(mfilename('fullpath'));
subfolder = ''; %include '/' at end of subfolder name
fname = 'HW6_pblm3';
simTime = '20';


%% System Definitions

% Ring Oscillator Parameters
% Based on the model w/ \dot{x} = -x + u  & y = alpha*tanh(beta*x)
tau = 1;
alpha = 1;
beta = mu/alpha;
x0 = 1;

tau1 = tau;
alpha1 = 2*alpha;
beta1 = beta/2;
x01 = x0;

tau2 = tau;
alpha2 = alpha;
beta2 = beta;
x02 = -2*x0;

tau3 = tau;
alpha3 = alpha;
beta3 = beta;
x03 = 3*x0;

% Coupling Matrix
K = toeplitz([0, -1, 0],[0, 0, -1]);



if openModel
    open(fname); % Don't need to open to run
end

% Auto Arrange
Simulink.BlockDiagram.arrangeSystem(gcs) %Auto Arrange
print(['-s', gcs], '-dpng',... % Print model to figure
    [cfolder, '\' subfolder, 'fig\', 'sim_model_', fname, '.png'])

if simulateModel
%% Simulate System
simConfig.SaveState = 'on';
simOut = sim(fname, 'SaveState', 'on', 'StartTime', '0', 'StopTime', simTime);

% Sim Data
Y_out = simOut.yout{1}.Values;
Xout = simOut.xout{1}.Values.Data; %Only works by grabbing states of first block (LTI_sys)
end

if plotResults
%% Plot Results

fig = figure;
hold on
plot(Y_out.Data(:,1))
plot(Y_out.Data(:,2))
plot(Y_out.Data(:,3))
title(['3-Stage Oscillator with \mu = ', string(mu)])
legend('y_1','y_2','y_3')

mu_str = string(1*(mu));
saveas(fig, string([cfolder, '\',subfolder, 'fig\', fname, '_results','_mu_',num2str(mu),'.png']))
% close all
end