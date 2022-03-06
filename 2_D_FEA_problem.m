clc; clear all;
%%
% An incomplete code for parallel plate capacitor is given below.
% A few major parameters are missing in the code. You will get hint of
% missing part once you run the program. Try to follow
% the steps taught in the class to compute and plot the value of
% voltage.
%% Equation to be solved: Poisson's equation
% form: ( Del(Del(phi)) = -rho/epsilon )
% Geometry of parallel plate capacitor is as follows:
% Voltage(y=1) = 10 volts
% ........................................................... y = 1
% . .
% . epsr=1 .
% . .
% ........................................................... y = 0
% x=0 Voltage(y=0) = 10 volts x=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%% Specifying limits of 2D geometry
Xlimits=[0 , 1];
Ylimits=[0 , 1];
sizeEx = 0.2; %size of element in x
sizeEy = 0.2; %size of element in y
%% Specify material properties
display('MATERIAL PARAMETERS');
BCvalue_top = input('voltage at upper plate= ');
BCvalue_bottom = input('voltage at upper plate= ');
epsilon_0 = 8.85e-12; % vaccum permittivity
epsilon_r = 1; % relative permittivity
epsilon = epsilon_0*epsilon_r; % permittivity of the dielectric
Rho= 0; % space charge in volume of dielectric
%% Meshing
% Assigning coordinates
[x,y] = meshgrid(Xlimits(1):sizeEx:Xlimits(2),Ylimits(1):sizeEy:Ylimits(2));
DT = delaunay(x,y); % Generate triangular mesh
M = size(DT,1); % number of elements in mesh
Nt = max(max(DT)); % total number of nodes in the mesh
nx = 1+ floor(Xlimits(2)- Xlimits(1))/sizeEx; % number of nodes in X direction
ny= 1+ floor( ( Ylimits(2)- Ylimits(1))/sizeEy ); %number of nodes in Y direction
% plot the mesh
triplot(DT,x,y); % Plotting mesh
axis equal tight
xlabel('x'); ylabel('y');
%% Mapping the X-Y coordinates of each local node in each element to correspinding global nodes
for p= 1:M
for i= 1:3
xe{p}(i) = x( DT(p,i)); % x coordinate of i'th local node in pth element
ye{p}(i) = y( DT(p,i)); % y coordinate of i'th local node in pth element
end
end
%%
% constants in shape function (p,q,r) -> S = p+qx+ry
% shape function is the parameter which is need to be defined for completion of the code
%% Computation of local force matrix
for p = 1:M
Area_tri = (1/2)*det([1, xe{p}(1), ye{p}(1) ;...
1, xe{p}(2), ye{p}(2) ;...
1, xe{p}(3), ye{p}(3) ]); % area of the triangular element
D = 2*Area_tri; % determinant
q1 = ye{p}(2) - ye{p}(3);
q2 = ye{p}(3) - ye{p}(1);
q3 = ye{p}(1) - ye{p}(2);
r1 = xe{p}(3) - xe{p}(2);
r2 = xe{p}(1) - xe{p}(3);
r3 = xe{p}(2) - xe{p}(1);
Ke_mat{p} = (epsilon/(2*D)) * [(q1*q1+r1*r1) (q1*q2+r1*r2) (q1*q3+r1*r3);
                               (q2*q1+r2*r1) (q2*q2+r2*r2) (q2*q3+r2*r3);
                               (q3*q1+r3*r1) (q3*q2+r3*r2) (q3*q3+r3*r3)];
Fe_mat{p} = -(Rho*D/6)*[1;1;1] ; % local force matrix for p'th element
end
%% Assembly of Global Stiffness Matrix
K = zeros(Nt,Nt); %Initializing global stiffness matrix
for p = 1:M %iterating over all elements
for i=1:3 %iterating over all the rows of local stiffness matrix
for j= 1:3 %iterating over all the columns of local stiffness matrix
K(DT(p,i), DT(p,j)) = K(DT(p,i), DT(p,j))+ Ke_mat{p}(i,j);
%Assembling the values at global node by mapping values at local nodes (thru all elements) into it
end
end
end
%% Assembly of Global force matrix
Force= zeros(Nt,1); %Initializing global force matrix
for p = 1:M
for i=1:3
Force(DT(p,i)) = Force(DT(p,i)) +Fe_mat{p}(i);
end
end
%% Applying Boundary conditions
BCnodes_top = []; % intializing boundary nodes at top of parallel plate
BCnodes_bottom = []; % intializing boundary nodes at bottom of parallel plate
for i= 1:Nt
if rem(i,ny)==0
BCnodes_top = [BCnodes_top, i]; % Boundary nodes at top electrode
end
if rem(i,ny)==1
BCnodes_bottom = [BCnodes_bottom, i]; % Boundary nodes at bottom electrode
end
end
for i = 1:length(BCnodes_top) % loop over all top boundary nodes
nBC = BCnodes_top(i); % boundary node
K(nBC,:) = 0; % make the n-th row of global K matrix zero
K(nBC,nBC) = 1; % make the diagonal term of the global K matrix 1
Force(nBC) = BCvalue_top; % assign the n-th row of Force to the known boundary value
end
for i = 1:length(BCnodes_bottom) % loop over all bottom boundary nodes
nBC = BCnodes_bottom(i); % boundary node
K(nBC,:) = 0; % make the n-th row of global K matrix zero
K(nBC,nBC) = 1; % make the diagonal term of the global K matrix 1
Force(nBC) = BCvalue_bottom; % assign the n-th row of Force to the known boundary value
end
%% Solution to FEM
%Equation is of form [K][Phi]= [F]
Phi = (inv(K))*Force; %nodal voltage values
Phiout = ([reshape(Phi,[nx,ny] )]); %rearanging nodal values in nx X ny matrix
%% Plotting solution
pcolor(x,y,Phiout)
shading interp
colorbar
colormap 'jet'
title('Voltage distribution')
xlabel('distance in x')
ylabel('distance in y')