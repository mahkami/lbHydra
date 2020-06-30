% Two Dimensional Example Postprocessor:
%   This postprocessor script:
%     1. Reads in the fluid densities from the given file
%     2. Extract a XY slice in the middle of Z axis as a sample
%     3. Calculate expected velocity distribution in this slice using analytical method
%     4. Graphically overlay the actual and analytical velocities for comparison
%   
%   The input file to analyze and the slice to be considered can be changed by 
%   altering this postprocessor.
%
%   To run this postprocessor and view the results, just type
%   "Postprocessor_TwoD" within Matlab.

%%%%%%%%%%%
%         %
%  Setup  %
%         %
%%%%%%%%%%%

clear all;

% Include Lattice-Boltzmann Matlab Scripts in search path
% Euler path
LBMatlab_Path = '/cluster/scratch/mahkami/Simulations/Permeability-real case/LBHMatlabToolbox';

% Mac path
LBMatlab_Path = '/Volumes/erdw_ifg_saar_home_01$/mahkami/Experimental Design/Simulations/lbHydra/Permeability/Overall Permeability/LBHMatlabToolbox';

add = pwd;

LBMatlab_Path = [ add '/LBHMatlabToolbox'];

addpath(LBMatlab_Path);

%%%%% Filenames %%%%%

% Simulation input files
coordFile = 'TwoDExample.nodeMap';

%  dat = [20000:20000:820000];
% dat = [2000:2000:30000 5000:5000:45000 55000:5000:80000 20000:20000:80000 120000:20000:180000 220000 250000:50000:550000];
dat = [820000];

for i = 1 : length(dat)
    m = int2str(dat(i));
    disp(['iteration' m])
    tic
    
    if ((dat(i) < 100000) && (dat(i)>= 10000))
        
        outputFile = ['out_0',m,'.lbData_bin'];
    elseif dat(i) < 10000
        outputFile = ['out_00',m,'.lbData_bin'];
    else
        
        outputFile = ['out_',m,'.lbData_bin'];
    end
    
   


% Output file to analyse
% outputFile = 'out_010000.lbData_bin';
% outputFile = 'out_800000.lbData_bin';
%%%%% Simulation Parameters %%%%%

% Load node coordinates
coords = load(coordFile);

%%%%%%%%%%%%%%%%%%%%%
%                   %
% Load output data  %
%                   %  
%%%%%%%%%%%%%%%%%%%%%

data = readLBData(outputFile);

% lbdata is set to zero in solid nodes.
indx = find(coords(:,5)==1);
data(indx,:) = 0;

%%
%%%%%%%%%%%%%%%%%%%%%%%
%                     %                                    
% Calculate velocity  %
%                     %  
%%%%%%%%%%%%%%%%%%%%%%%
c = d3q19LatticeDirections; % get the d2q9 lattice directions
v = calculateVelocities(data,c);

% lattice dimensions
x_node = 1 + max(coords(:,1));
y_node = 1 + max(coords(:,2));
z_node = 1 + max(coords(:,3));

% Velocity components
v_x = v(:,1);
max(max(v_x(:,:)))
v_y = v(:,2);
v_z = v(:,3);

% reshaping the velocity components
v_x = reshape(v_x,x_node,y_node,z_node);
v_y = reshape(v_y,x_node,y_node,z_node);
v_z = reshape(v_z,x_node,y_node,z_node);


%%%%%%%%%%%%%%%%%%%%%%%
%                     %                                    
% Calculate Pressure  %
%                     %  
%%%%%%%%%%%%%%%%%%%%%%%

rho = calculateDensities(data);
rho = reshape(rho,x_node,y_node,z_node);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %                                    
% Calculate Permeability  %
%                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k = u*mu*dx/dp

kinematicViscosity = 0.16666667; % kinematic viscosity in lattice units. 

% inp_k = zeros(y_node,4);
% inp.u_in1  = v_x(:,:,1);  % inlet velocity
% inp.u_in2  = v_x(:,:,2);  % inlet velocity
% inp.u_out1 = v_x(end,:,1);  % outlet velocity
% inp.u_out2 = v_x(end,:,2);  % outlet velocity
% inp.p_in1  = rho(1,:,1); % inlet pressure
% inp.p_in2  = rho(1,:,2); % inlet pressure
% inp.p_out1 = rho(end,:,1); % outlet pressure
% inp.p_out2 = rho(end,:,2); % outlet pressure

% inp.u_in = mean(inp.u_in1(inp.u_in1>0))+mean(inp.u_in2(inp.u_in>0));  % average inlet velocity
% inp.u_out = mean(inp.u_out1(inp.u_out1>0))+mean(inp.u_out2(inp.u_out>0));  % average outlet velocity
% 
% inp.p_in = mean(inp.p_in1(inp.u_in1>0))+mean(inp.p_in2(inp.p_in>0));  % average inlet pressure
% inp.p_out = mean(inp.p_out1(inp.u_out1>0))+mean(inp.p_out2(inp.p_out>0));  % average outlet pressure


inp.u_in  = v_x(:,:,1:2);  % inlet velocity

inp.p_in1  = rho(1,:,1:2); % inlet pressure

inp.p_out1 = rho(end,:,1:2); % outlet pressure


inp.u = mean2(inp.u_in(inp.u_in>0));  % average total velocity
inp.p_in = mean2(inp.p_in1(inp.p_in1>0)); % average inlet pressure
inp.p_out = mean2(inp.p_out1(inp.p_out1>0)); % average inlet pressure

k(i) = -inp.u * kinematicViscosity * x_node/(inp.p_out-inp.p_in);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         %                                    
% Unit Conversion         %
%                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cell_size = [70 0.5]; % Cell width and depth in mm
% % Calculating Reynolds number
% % Re = u*Dh/kinemaatic viscosity, Dh = Cross section are/Perimiter
% 
% DH.LB = 4*y_node/(2*(4+y_node));
% Re = inp.u * DH.LB/kinematicViscosity;
% 
% % u/u_lb = v/v_lb * DH_lb/DH
% DH.exp = 1e-3 * cell_size(1)*cell_size(2)/2*(cell_size(1)+cell_size(2));
% water_viscosity = 1e-6;
% 
% u_conv = water_viscosity/kinematicViscosity * DH.LB/DH.exp;
% u_conv = u_conv * 1e3 * 60; 
% 
% v_x = v_x * u_conv;
% v_y = v_y * u_conv;

% Normalizing velocity with injection velocity
% load('inection_rate');
% v_x = v_x /inj_vel;
% v_y = v_y /inj_vel;

% Normalizing with the velocity in the inlet
% v_x = v_x /max(v_x(2,:,2));
% v_y = v_y /max(v_y(2,:,2));


% Normalizing with maximum velocity
% 
% v_x = v_x /max(v_x(2,:,2));
% v_y = v_y /max(v_y(2,:,2));


%% Calculating permeability 
load('viscosity.mat');
load('deltarho.mat');
% average velocity 
av_vel = mean(v(10:end-9,1));

perm_lb =  -kinematicViscosity*av_vel/(deltaRho/x_node) % permeability in lbm unit

% % Permeability conversion
% one lattice is 0.05 mm >>>> one lattice is 5e-5 m
% (lb?2) * (m/lb)?2 = m?2
conv = 5e-5;
perm =  perm_lb* conv^2 
figure(120);
plot(dat(i),perm, '*','color','red');
% saveas(gcf,'Permeability','epsc')
savefig(gcf,'Permeability')
hold on

perm_out(i) = perm; 
iteration(i) = dat(i);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           %                                    
% Plot and compare results  %
%                           %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  load('x-dimension.mat')
%  load('y-dimension.mat')
%  
 tt = num2str(i);
% vtkwrite(['vtk' tt'.vtk'], 'structured_grid',x_node,y_node,z_node,'vectors','velocity',v_x,v_y,v_z)
% vtkwrite(['vtk' m '.vtk'], 'structured_points', 'velocity', v_x(:,:,2:7))
% Mat2VTK(['Velocity' m '.vtk'] ,abs((v_x(:,:,2:7))),'binary')
% Mat2VTK(['Pressure' m '.vtk'] ,abs((rho(:,:,2:7))),'binary')
% animation(v_x,length(dat),i);

%y = length(v_x(1,:,2));

%for jj = 1:10
 %   hfig = figure('visible','off');
 %   peak = findpeaks(v_x(1+100*(jj-1),:,2));
 %   x = linspace(0,y,length(peak));
 %   plot(x,peak)
 %   st = num2str(jj*100);
 %   legend(st)
 %   axis([0 y 0 100])
 %   M(jj) = getframe(hfig);
 %   close
%end

%figure(12), subplot(1,length(dat),i) 
%axes('position',[0 0 1 1],'visible','off');
% movie(M,2,1)

 


% subplot(1,length(dat),i)
% figure()
% imagesc(abs(v_y(:,:,4)))
% title(['Velocity in y direction' m],'fontsize',18)
% colormap(jet);
% colorbar;
% caxis([0 5e-8])
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% % saveas(gcf,'Velocity in y direction','epsc')
% savefig(gcf,['Velocity in y direction' m])
% hold on
% 

% % subplot(1,length(dat),i)
% figure()
% imagesc(abs(v_x(:,:,4)))
% title(['Velocity in flow direction' m],'fontsize',18)
% colormap(jet);
% colorbar;
% % caxis([0 max(max(v_x(:,:,2)))])
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% saveas(gcf,'Velocity in flow direction','epsc')
% savefig(gcf,['Velocity in flow direction' m])
% % hold on


% % subplot(1,length(dat),i)
% figure()
% imagesc(abs(v_z(:,:,4)))
% title(['Velocity in z direction' m],'fontsize',18)
% colormap(jet);
% colorbar;
% % caxis([0 max(max(v_x(:,:,2)))])
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% % saveas(gcf,'Velocity in flow direction','epsc')
% savefig(gcf,['Velocity in z direction' m])
% % hold on



% 
% % figure(3)
% % imagesc(v_x(2:300,:,2))
% % subplot(1,length(dat),i)
% % title(['Velocity in x direction - Near Injection point',' Out',m],'fontsize',18)
% % colormap(jet);
% % colorbar;
% % % caxis([0 0.003])
% % % caxis([0 max(max(v_x(2:100,:)))])
% % xlabel('Cell width','fontsize',18)
% % ylabel('Cell length','fontsize',18)
% % saveas(gcf,'Velocity in x direction - Near Injection point','epsc')
% % savefig(gcf,'Velocity in x direction - Near Injection point')
% 
% 
% figure()
% % subplot(1,length(dat),i)
% quiver(v_y(1:5:end,1:5:end,4),v_x(1:5:end,1:5:end,4),5);
% % quiver(v_y(:,:,2),v_x(:,:,2));
% title(['velocity vectors' m],'fontsize',18)
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% % set(gca,'Ydir','reverse')
% % saveas(gcf,'velocity vectors','epsc')
% savefig(gcf,['velocity vectors' m])
% % hold on
% 
% 
% 

% % subplot(1,length(dat),i)
% figure()
% imagesc(rho(:,:,4))
% title(['Densities' m],'fontsize',18)
% colormap(jet);
% colorbar;
% % caxis([0 0.003])
% % caxis([1 max(max(rho(:,:,2)))])
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% % saveas(gcf,'Densities')
% savefig(gcf,['Densities' m])
% % hold on
% 
% figure(6)
% subplot(1,length(dat),i)
% plot(v_x(round(x_node/10),:,2))
% title(['velocity profile 10%',m],'fontsize',18)
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% saveas(gcf,'velocity profile 10%')
% savefig(gcf,'velocity profile 10%')
% hold on
% 
% figure(7)
% subplot(1,length(dat),i)
% plot(v_x(round(9*x_node/10),:,2))
% title(['velocity profile 90%',m],'fontsize',18)
% xlabel('Cell width','fontsize',18)
% ylabel('Cell length','fontsize',18)
% saveas(gcf,'velocity profile 90%')
% savefig(gcf,'velocity profile 90%')
% hold on
% 
% figure(8)
% % subplot(1,length(dat),i)
% plot(dat(i),k(i),'*','color','r')
% title(['Permeability'],'fontsize',18)
% xlabel('Iteration','fontsize',18)
% ylabel('Permeability [mDa]','fontsize',18)
% saveas(gcf,'Permeability')
% savefig(gcf,'Permeability')
% hold on
toc


end
save('Permeability','perm_out')
save('iteration','iteration')

disp('PostProcessor finished')
