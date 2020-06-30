function dom_4x = Lower_Part_Final()

% The script designs a system of porous fracture

% the inputs are:
% Minimum pixels per smallest capillary
% smallest capillary size
% camera resolutions

% The target of this cell is to observe:
% high permeable to low permable porous flow
% Porous to dead-end fracture flow
% Porous to staright-run fracture flow
% Fracture (dead-end) to Fracture(straigt-run) flow





clear all


bounce    = 1; % will map 1 to LBD2Q9 bounceback model
interior  = 0; % 2 to LBD2Q9 model
low_dens  = 3; % mapping 3 to density on lower boundary
high_dens = 4; % mapping 4 to density on upper boundary


% cell_l = round((y_res*min_app)/pixPap); % cell width
% cell_w = round(cell_l*x_res/y_res);     % cell length

% Manual input of cell size to match the camera. We place 1 cell for each
% 0.1 mm.

cell_l = 84;
cell_w = 100;

% Defining the geometry with dimensions of  w   l
% Each cell has the dimension of 0.1*0.1 mm?2

% x_nodes = cell_w/min_app;
% y_nodes = cell_l/min_app;

x_nodes = cell_w*10;
y_nodes = cell_l*10;

dom_p1 = interior*ones(x_nodes,y_nodes);
dom_p2 = interior*ones(x_nodes,y_nodes);
dom    = interior*ones(x_nodes,y_nodes);

% 0.2 mm porous media - pin size 0.8 mm
t = 0; % counter

pins = 8; %0.8 mm size of the cell
space_1 = 2; % 0.2 mm spacing between pins
space_2 = 3; % 0.3 mm spacing between pins

i = 1;
while i <= x_nodes-1
    for sp = 1 : space_1
        if i ~= x_nodes
            i = i+1;
        end
        
        
    end
    sp = 0;
    
    for pin = 1 : pins
        
        if i ~= x_nodes
            t = t+1;
            x_p1(t) = i;
            i = i+1;
        end
        
    end
    pin = 0;
end

t = 0;
i = 1;

while i <= y_nodes-1
    for sp = 1 : space_1
        if i ~= y_nodes
            i = i+1;
        end
    end
    sp = 0;
    
    for pin = 1 : pins
        
        if i ~= y_nodes
            t = t+1;
            y_p1(t) = i;
            i = i+1;
        end
        
    end
    pin = 0;
end

dom_p1(x_p1,y_p1) = bounce;

% 0.3 mm porous media - pin size 0.8 mm

t= 0;  % counter
i = 1;
while i <= x_nodes-1
    for sp = 1 : space_2
        if i ~= (x_nodes-1)
            i = i+1;
        end
        
    end
    sp = 0;
    
    for pin = 1 : pins
        
        if i ~= x_nodes
            t = t+1;
            x_p2(t) = i;
            i = i+1;
        end
        
    end
    pin = 0;
end

t = 0;
i = 1;
while i <= y_nodes-1
    for sp = 1 : space_2
        if i ~= y_nodes
            i = i+1;
        end
        
    end
    sp = 0;
    
    for pin = 1 : pins
        
        
        if i ~= y_nodes
            t = t+1;
            y_p2(t) = i;
            i = i+1;
        end
        
    end
    pin = 0;
end

dom_p2(x_p2,y_p2) = bounce;


% Putting two porous media in the cell

dom(:,1:round(y_nodes)/2+10) = dom_p1(:,1:round(y_nodes)/2+10);
dom(:,round(y_nodes)/2+11:y_nodes) = dom_p2(:,round(y_nodes)/2+11:y_nodes);


% introducing the fractures
[x y] = size(dom);


% staright through flow fractures

% fracture in zone 1
frac_high = 40;

dom(:, round(y/2 - y/7 -20):round(y/2 - y/7))= interior;

% fracture in zone 2

dom(:, round(y/2 +y/7):round(y/2 + y/7 + 20))= interior;


[x y] = size(dom);

% oriented isolated fracture
% each point is 0.1 mm

% fracture in zone 1


t = 0;
jj = y/2;
while jj > pins+20
    jj = 152-t*(pins+space_1);
    
    frac_y = (jj - 20:jj);
    frac_x= (t*(6*pins+6*space_1)+1 :(t+1)*(6*pins+6*space_1)+1);
    t = t+1;
    dom(frac_x,frac_y) = interior;
       
end

jj = 0;
t=0;
while jj < y-pins-20
    jj = 686 +t*(pins+space_2);
    
    frac_y = (jj:jj+20+1);
    frac_x= (t*(6*pins+6*space_2)+1 :(t+1)*(6*pins+6*space_2)+1);
    t = t+1;
    dom(frac_x,frac_y) = interior;
       
end

% adding 10% of additional length for uniform flow in each side

% 5% for the oriented section and 5% void space

temp = interior*zeros(x + round(x/10),y);

temp(round(x/20)+1:end-round(x/20),:) = dom;

% calculating the slope of the oriented section

% diameter of input 11.4 mm = 114 pixel
% slope is the ratio of x to y so we can know the changing ratio of y cells
% with the change in 1 x cell

% x change in x direction is 20% of x length and y direction is y -
% injection length

% slope = round((x/6)\(y/2-57));
%
%
% if slope > 0

dom = temp;

%Updating the geometry

[x y] = size(dom);

inj = 40;
% adding 4 mm for wall thickness on left and right side

temp = imresize(logical(dom),2,'box');
temp = double(temp);
dom = temp;

dom(1,:) = bounce;
dom(end,:) = bounce;

dom = flip(dom,1);

[x y] = size(dom);

dom(1  ,y/2-inj:y/2+inj) = low_dens;
dom(end,y/2-inj:y/2+inj) = high_dens;

[x y] = size(dom);
% Updating the cell size in mm
cell_l = cell_l;  % 4 mm in each side

% 50% added for void space, 30mm added for threads
cell_w = cell_w + 0.1*cell_w;

%%
x_disc = linspace(0,cell_w,x);
y_disc = linspace(0,cell_l,y);

% saving dimensions

save('x-dimension', 'x_disc');
save('y-dimension', 'y_disc');

%% Plotting the geometry

figure, imagesc(y_disc,x_disc,dom);
title('Cell Geometry','fontsize',18)
% colormap(jet);
xlabel('Cell width (mm)','fontsize',18)
ylabel('Cell length (mm)','fontsize',18)
saveas(gcf,'Cell Geometry','epsc')
savefig(gcf,'Cell Geometry')
disp('Geometry finished')

%% Converting the Geometry to 3D
cell_h = 0.3; % 4mm cell height
nz = 8; % 6 nodes in z direction

z_disc = linspace(0,cell_h,nz);
save('z-dimension', 'z_disc');

[x y] = size(dom)
nz
temp = zeros(x,y,nz);

interior  = 2;
dom(dom<1) = interior;


for ii = 2 : (nz -1)
    temp(:,:,ii) = dom;
end

temp(:,:,1) = bounce;
temp(:,:,end) = bounce;


dom_4x = temp;




disp('Geometry finished')

end

