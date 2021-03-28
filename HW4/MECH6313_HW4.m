% MECH 6313 - Homework 4
% Jonas Wagner
% 2021-03-26

clear
close all

%% Simulink Creation In Code
% Specify the name of the model to create
fname = 'MECH6313_HW4_pblm3_model';

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


%% Problem 3
% Saturation Block
a = -1;
b = 1;

% Linearized System
A = [0, 1;
     0, 0];
B = [0;
     1];
C = eye(2);
D = [0;
     0];
lti_sys = ss(A,B,C,D);
x0 = [1;
     -1];


% Controller Gain
K = place(A,B, [-1+j,-1-j]);
A_BK = A - B * K;
eig_A_BK = eig(A_BK);


% Creat Simulink Model
new_system(fname);

% Create Simiple Input
add_block('simulink/Sources/In1', [gcs, '/In']);

% Create Sum block
add_block('simulink/Commonly Used Blocks/Sum', [gcs, '/Sum'],...
    'inputs', '|+-');
add_line(fname, 'In/1', 'Sum/1','autorouting','on');

% Saturation Block
add_block('simulink/Commonly Used Blocks/Saturation', [gcs, '/Saturation'], ...
    'LowerLimit','a', ...
    'UpperLimit','b');
add_line(fname, 'Sum/1', 'Saturation/1');

% State-Space System
add_block('cstblocks/LTI System', [gcs, '/LTI_sys'],...
    'sys','lti_sys',...
    'IC', 'x0');
add_line(fname, 'Saturation/1', 'LTI_sys/1');

% Controller
add_block('simulink/Commonly Used Blocks/Gain', [gcs, '/Controller'],...
    'Gain', 'K',...
    'Multiplication', 'Matrix(K*u)');
add_line(fname, 'LTI_sys/1', 'Controller/1');
add_line(fname, 'Controller/1', 'Sum/2');

% Create Simple Scope/Output
add_block('simulink/Sinks/Out1', [gcs, '/Out']);
add_line(fname, 'LTI_sys/1', 'Out/1');


% Auto Arrange
Simulink.BlockDiagram.arrangeSystem(fname) %Auto Arrange

% Save System
save_system(fname);
open(fname);

% Simulate System
simConfig.SaveState = 'on';
simOut = sim(fname, simConfig);

% Sim Data
Xout = simOut.xout{1}.Values.Data;

% Plot Outputs
figure
plot(Xout(:,1))
hold on
plot(Xout(:,2))

