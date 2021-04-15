function param_est(block)
% Level-2 MATLAB file S-Function for continuous time variable step demo.

%   Copyright 1990-2009 The MathWorks, Inc.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of input and output ports
  % Y
  % Psi
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = [3, 3];
  block.InputPort(1).DirectFeedthrough = true;
  block.InputPort(2).Dimensions        = [1, 3];
  block.InputPort(2).DirectFeedthrough = true;
  
  
  block.OutputPort(1).Dimensions       = 1;
  
  %% Set block sample time to variable sample time
  block.SampleTimes = [-2 0];
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions); 
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Update',                  @Update);  

%endfunction

function DoPostPropSetup(block)

  %% Setup Dwork
  block.NumDworks = 1;
  block.Dwork(1).Name = 'X'; 
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;

%endfunction

function InitConditions(block)

  %% Initialize Dwork
  block.Dwork(1).Data = 0;
  
%endfunction

function Output(block)

  block.OutputPort(1).Data = Theta_Dot;
%   block.Dwork(1).Data;
  
  %% Set the next hit for this block 
  block.NextTimeHit = block.CurrentTime + block.InputPort(1).Data(2);
  
%endfunction

function Update(block)
  
  Psi = block.InputPort(1).Data;
  Y = block.InputPort(2).Data;
  Theta_dot = Psi' * Psi - Psi * Y;
  block.Dwork(1).Data = Theta_dot;
  
  
%endfunction

