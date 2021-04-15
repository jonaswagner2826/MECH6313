% MECH 6313 - HW5
% Jonas Wagner
% 2021-04-08
% 

% Only Problem 1b
clear
close all

%% Settings
generateModel = true;
openModel = false;
simulateModel = true;
plotResults = true;

% Name of the simulink model
[cfolder,~,~] = fileparts(mfilename('fullpath'));
subfolder = ''; %include / at end of subfolder
fname = 'pblm1b';


%% System Definitions

% Physical Parameters
m = 20;
b = 0.1;
k = 5;

% Matrix Definitions
A = [-b/m, -k/m;
     1, 0];
B = [1/m;
     0];
C = [0, 1];
D = 0;
lti_sys = ss(A,B,C,D);
x0 = [1;
     -1];

% Filter TF
lambda_0 = 20;
lambda_1 = 2 * 1 * lambda_0;
filter_1 = tf([1], [1, lambda_1, lambda_0]);
% filter_1 = [filter_1; filter_1];
filter_2 = tf([1, 0], [1, lambda_1, lambda_0]);
% filter_2 = [filter_2; filter_2];
filter_3 = - filter_1;


if generateModel
%% Simulink Creation In Code
% Simulink Settings ----------------------
% Get the current configuration
cfg = Simulink.fileGenControl('getConfig');
% Changes Code Save Location
cfg.CacheFolder = [pwd, '\', subfolder];
cfg.CodeGenFolder = [pwd, '\', subfolder];
cfg.CodeGenFolderStructure = 'TargetEnvironmentSubfolder';
% Apply new Config
Simulink.fileGenControl('setConfig', 'config', cfg, 'keepPreviousPath',true, 'createDir',true);

% Check if the file already exists and delete it if it does
if exist(fname,'file') == 4
    % If it does then check whether it's open
    if bdIsLoaded(fname)
        % If it is then close it (without saving!)
        close_system(fname,0)
    end
    % delete the file
    delete([fname,'.slx']);
end

% Create Simulink Model
new_system;

% Create Simiple Input
add_block('simulink/Sources/Sine Wave', [gcs, '/In']);

% Plant
add_block('cstblocks/LTI System', [gcs, '/LTI_sys'],...
    'sys','lti_sys',...
    'IC', 'x0');
u = add_line(gcs, 'In/1', 'LTI_sys/1');
set_param(u, 'Name', 'u')


% Filter 1
add_block('cstblocks/LTI System', [gcs, '/filter_1'],...
    'sys','filter_1');
add_line(gcs, 'In/1', 'filter_1/1');


% Filter 2
add_block('cstblocks/LTI System', [gcs, '/filter_2'],...
    'sys','filter_2');
add_line(gcs, 'LTI_sys/1', 'filter_2/1');


% Filter 1
add_block('cstblocks/LTI System', [gcs, '/filter_3'],...
    'sys','filter_3');
add_line(gcs, 'LTI_sys/1', 'filter_3/1');


% mux
add_block('simulink/Commonly Used Blocks/Mux', [gcs,'/mux'],...
    'inputs', '3');
psi_1 = add_line(gcs, 'filter_1/1', 'mux/1');
psi_2 = add_line(gcs, 'filter_2/1', 'mux/2');
psi_3 = add_line(gcs, 'filter_3/1', 'mux/3');
set_param(psi_1, 'Name', 'psi_1')
set_param(psi_2, 'Name', 'psi_2')
set_param(psi_3, 'Name', 'psi_3')


% Parameter Estimator
% -------------------------
% Negative Gain Block
add_block('simulink/Commonly Used Blocks/Gain', [gcs,'/Psi_neg'],...
    'gain', '-1');
Psi = add_line(gcs, 'mux/1', 'Psi_neg/1');
set_param(Psi, 'Name', 'Psi')



% Product Block
add_block('simulink/Math Operations/Product', [gcs, '/Mult1'],...
    'inputs', '**',...
    'Multiplication', 'Matrix(*)');
Psi_neg = add_line(gcs, 'Psi_neg/1', 'Mult1/1');
set_param(Psi_neg, 'Name', 'Psi_neg')
Y = add_line(gcs, 'LTI_sys/1', 'Mult1/2');
set_param(Y, 'Name', 'Y')


% Create Sum block
add_block('simulink/Commonly Used Blocks/Sum', [gcs, '/Sum1'],...
    'inputs', '|+-');
theta_m = add_line(gcs, 'Mult1/1', 'Sum1/1');
set_param(theta_m, 'Name', 'theta_m')

% Creates Integrator
add_block('simulink/Commonly Used Blocks/Integrator', [gcs, '/Int1']);
theta_error = add_line(gcs, 'Sum1/1', 'Int1/1');
set_param(theta_error, 'Name', 'theta_error');

% Create Transpose
add_block('simulink/Math Operations/Math Function', [gcs, '/Trans1'],...
    'Function', 'transpose')
add_line(gcs, 'mux/1', 'Trans1/1')

% Product Block (feedback)
add_block('simulink/Math Operations/Product', [gcs, '/Mult2'],...
    'inputs', '***',...
    'Multiplication', 'Matrix(*)');
Psi = add_line(gcs, 'mux/1', 'Mult2/1');
set_param(Psi, 'Name', 'Psi')
add_line(gcs, 'Trans1/1', 'Mult2/2');
% add_param(Psi_T, 'Name', 'Psi_T')
theta_est = add_line(gcs, 'Int1/1', 'Mult2/3');
set_param(theta_est, 'Name', 'theta_est')

add_line(gcs, 'Mult2/1', 'Sum1/2')


% Create Simple Scope/Output
add_block('simulink/Sinks/Out1', [gcs, '/Y']);
y = add_line(gcs, 'LTI_sys/1', 'Y/1');
set_param(y, 'Name', 'y');
add_block('simulink/Sinks/Out1', [gcs, '/Psi']);
Psi = add_line(gcs, 'mux/1', 'Psi/1');
set_param(Psi, 'Name', 'Psi');
add_block('simulink/Sinks/Out1', [gcs, '/Theta_est']);
Theta_est = add_line(gcs, 'Int1/1', 'Theta_est/1');
set_param(Theta_est, 'Name', 'Theta_est');



% https://www.mathworks.com/help/simulink/ug/creating-an-example-model-that-uses-a-matlab-function-block.html
% S-block implimentation... ran out of time for implimenatation
% add_block('simulink/User-Defined Functions/Level-2 MATLAB S-Function', [gcs, '/Est'],...
%     'FunctionName', 'param_est')
% add_line(gcs, 'Mult1/1', 'Est/1');
% add_line(gcs, 'LTI_sys/1', 'Est/2');


% Auto Arrange
Simulink.BlockDiagram.arrangeSystem(gcs) %Auto Arrange

%% Save and Open System
save_system(gcs,[cfolder, '\', subfolder, fname]);
print(['-s', gcs], '-dpng',... % Print model to figure
    [cfolder, '\' subfolder, 'fig\', 'sim_model_', fname, '.png'])

end
if openModel
    open(fname); % Don't need to open to run
end

if simulateModel
%% Simulate System
simConfig.SaveState = 'on';
simOut = sim(fname, 'SaveState', 'on', 'StartTime', '0', 'StopTime', '50');

% Sim Data
Y_out = simOut.yout{1}.Values;
Psi_out = simOut.yout{2}.Values;
Theta_out = simOut.yout{3}.Values;
Xout = simOut.xout{1}.Values.Data; %Only works by grabbing states of first block (LTI_sys)
end

if plotResults
%% Plot Results

fig = figure;
subplot(3,1,1)
plot(Y_out.Data(:,1))
title('Y')


subplot(3,1,2)
plot(Psi_out.Data(:,1))
hold on
plot(Psi_out.Data(:,2))
hold on
plot(Psi_out.Data(:,3))
title('Psi')

subplot(3,1,3)
plot(Theta_out.Data(:,1))
hold on
plot(Theta_out.Data(:,2))
hold on
plot(Theta_out.Data(:,3))
title('Theta_est')

% fig = figure;
% plot(Xout(:,1))
% hold on
% plot(Xout(:,2))
% legend('X_1', 'X_2')
% title('Saturated Double Integration Response while stabalized')
saveas(fig, [cfolder, '\',subfolder, 'fig\', fname,'_results.png'])
close all
end
