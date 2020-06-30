
%  Two Dimensional Example Preprocessor:
%    This preprocessor script:
%       1. Defines the material models for each lattice node.
%       2. Creates the lattice graph and maps each node to its neighbors
%       3. Generates the initial conditions for the lattice, a uniform density gradient 
%          (from 0.9975 to 1.0025) along the z axis.
%       4. Outputs this information to the appropriate files, which are recorded in a
%          configuration file, "TwoD.simulationInput".
% 
%    The input files generated can be altered with a text editor or by editing 
%    this preprocessor script. 
%
%    To run the simulation and view the results
%       1. Copy the simulation binary "generalLBSimulation" into the current directory.
%       2. Start Matlab. 
%       3. Run this preprocessor script:
%            Type "Preprocessor_TwoD" within Matlab and hit return. 
%       4. Run the lattice-Boltzmann simulation:
%            Either type "!./generalLBSimulation < TwoDExample.simulationInput" within Matlab, or 
%            "./generalLBSimulation < TwoDExample.simulationInput" from the command line. 
%       5. The program runs for the specified maximum number of timesteps (10,000) and produces a series 
%          of output files, eg. "out_010000.lbData", at regular intervals as specified by the output 
%          frequency variable (every 1,000 timesteps).
%       6. After the program has finished, run the postprocessor script by typing 
%          "Postprocessor_TwoD" within Matlab.  

%%%%%%%%%%%
%         %                                    
%  Setup  %
%         %  
%%%%%%%%%%%

clear all;

% Add Lattice-Boltzmann Toolbox to Matlab search path
% Euler path
% LBMatlab_Path = '/cluster/scratch/mahkami/Simulations/Permeability-real case/LBHMatlabToolbox';
% mac path
% LBMatlab_Path = '/Users/mahkami/Desktop/Scratch/Simulations/Permeability-real case/LBHMatlabToolbox';
addpath(LBMatlab_Path);

%%%%% Filenames %%%%%

% Input files
inputConfigurationFile = 'TwoDExample.simulationInput';
nodeMapFile = 'TwoDExample.nodeMap';
neighborMapFile = 'TwoDExample.neighborMap';
matModelDataFile = 'TwoDExample.matModelData';
initialLBDataFile= 'TwoDExample.lbData_bin'; 

% Output file prefix eg. "out_010000.lbData" (for a single phase)
outputFilePrefix = 'out';

%%%%% Simulation Parameters %%%%%

% Lattice Type
latticeType = 'D3Q19';

% Lattice dimensions
% nX = 20; % size of X-dimension (rows in matlab)
% nY = 80; % size of Y-dimension (columns in matlab)
% nZ = 1; % size of Z-dimension

% Physical parameters 
deltaRho = 0.00005; % Initial density difference, average density = 1.0
kinematicViscosity = 0.16666667; % kinematic viscosity in lattice units. 
Collision_frequency = 1/(3*kinematicViscosity+0.5);

save('viscosity','kinematicViscosity');
save('deltarho','deltaRho');

% Duration and output
maxTimeSteps = 200000;
outputFrequency = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       %                                    
% Stage  One  - Define Material Models  %
%                                       %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Defining material models')

%%%%%% Define material models %%%%%%
% Model 1: Bounceback boundary model on edges
materialModels(1).name = 'd3q19BouncebackModel';   

% Model 2: d3q19Model { kinematic_viscosity,collision_frequency } in
% interior
materialModels(2).name = 'd3q19Model';  
materialModels(2).params.kinematic_viscosity = kinematicViscosity;
materialModels(2).params.collision_frequency = Collision_frequency;

% % Model 3: d3q19ImmiscibleVelocityControlledModel { boundary_norm,velocity } on upper boundary
% 
% materialModels(3).name = 'd3q19ImmiscibleVelocityControlledModel'; 
% materialModels(3).params.velocity = [0.000001 0 0];
% materialModels(3).params.boundary_norm = '-1 0.0 0.0';

% inj_vel = materialModels(3).params.velocity(1);

% save('inection_rate','inj_vel')

% % Model 3: d3q19DensityControlledModel { boundary_norm,boundary_dir,kinematic_viscosity,rho } on lower boundary
% High pressure boundary
materialModels(3).name = 'd3q19DensityControlledModel';  
materialModels(3).params.kinematic_viscosity = kinematicViscosity;
materialModels(3).params.rho = 1; % - 0.5*deltaRho;
materialModels(3).params.boundary_norm = '1.0 0.0 0.0';
% 
% % % Model 4:  d3q19DensityControlledModel { boundary_norm,boundary_dir,kinematic_viscosity,rho } on upper boundary
% low pressure boundary
 materialModels(4).name = 'd3q19DensityControlledModel'; 
 materialModels(4).params.kinematic_viscosity = kinematicViscosity;
 materialModels(4).params.rho = 1 - deltaRho;
 materialModels(4).params.boundary_norm = '-1.0 0.0 0.0';

% Model 4: d3q19ImmiscibleVelocityControlledModel { boundary_norm,velocity } on upper boundary

% materialModels(4).name = 'd3q19ImmiscibleVelocityControlledModel'; 
% materialModels(4).params.velocity = [0.000022 0 0];
% materialModels(4).params.boundary_norm = '-1 0.0 0.0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %             
% Stage  Two  -  Generate Lattice / Assign models to nodes  %
%                                                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Assigning material models to nodes')

%%%%%% Assign material models to nodes %%%%%%
% 0 indicates a hole (absence of a node) and values greater than 0 are mapped to the
% material models given above.

% node_input = 'voxels';
% load('voxels.mat');

% This is the place to load the geometry, it can be manualy produced here 
% or it can be loaded from an external file orfrom another matlab script
voxels = Lower_Part_Final();


% a simple cube gemetry with size of 100*100*100, low in x direction

% voxels = 2*ones(100,100,100);
% voxels(:,1,:) = 1;
% voxels(:,end,:) = 1;
% voxels(:,:,1) = 1;
% voxels(:,:,end) = 1;
% voxels(1,1:end-1,1:end-1) = 3; 
% voxels(end,1:end-1,1:end-1) = 4; 

% voxels = voxels(1:1600,:,:);

% voxels = ones(nX,nY); % will map 1 to LBD2Q9 bounceback model
% voxels(2:end-1,:) = 2; % mapping 2 to LBD2Q9 model
% voxels(2:end-1,1) = 3; % mapping 3 to density on lower boundary
% voxels(2:end-1,end) = 4; % mapping 4 to density on upper boundary

%%%%%% Determine node connectivity and make node-model map %%%%%%
c =  d3q19LatticeDirections; % get the d2q9 lattice directions
[~,nodeMap] = generateNodeNeighborList(voxels,c);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Stage Three  -  Generate Initial Conditions  %
%                                              %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Generating initial conditions')

%%%%%% Create density (pressure) gradient along y axis %%%%%
minX = min(nodeMap(:,1)); maxX =  max(nodeMap(:,1));
densities(nodeMap(:,4)) = 1.0 - deltaRho*(minX-nodeMap(:,1))/(minX-maxX) ;

%%%%%% Convert densities to lattice-Boltzmann data %%%%%
% F is an Nx9 matrix describing the initial state of the lattice Boltzmman lattice.
disp('calculate densities')

w = d3q19LatticeWeights;    % get the d2q9 lattice weights
F = calculateLBData(densities,[],w); % generate data for d2q9 lattice


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %                                    
% Stage Four   -  Write Data Files  %
%                                   %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Writing data files')

% dlmwrite(nodeMapFile, nodeMap, 'precision', 9);         % Write node -> coordinate map
% dlmwrite(neighborMapFile, neighborMap, 'precision', 9); % Write node -> list of neighbors map

saveascii(nodeMap,nodeMapFile,{','},{'%d'});
% saveascii(neighborMap,neighborMapFile,{','},{'%d'});

writeLBData(initialLBDataFile,F);                       % Write initial fluid density distribution

writeMaterialModelFile(matModelDataFile,materialModels); % Write material model data

%%%%%% Write simulation input file %%%%%% 
% set input parameters
simulationInput.lattice_type = latticeType;
simulationInput.output_file_prefix = outputFilePrefix;
simulationInput.initial_lbdata_file = initialLBDataFile;
simulationInput.max_timesteps = maxTimeSteps;
simulationInput.output_frequency = outputFrequency;
simulationInput.node_connectivity_map_file = neighborMapFile;
simulationInput.node_to_material_model_map_file = nodeMapFile;
simulationInput.material_model_data_file = matModelDataFile;

writeSimulationInputFile(inputConfigurationFile, simulationInput); % Write simulation input file

%%%%%%%%%%
%        %                                    
%  END   %
%        %  
%%%%%%%%%%

disp('Finished')

%%%%% Done %%%%%

