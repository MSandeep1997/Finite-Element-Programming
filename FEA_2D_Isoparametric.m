%This is a MATLAB programme to solve the 2-D FEM programme using
%Isoparametric elements
%NAME: Sandeep Mohan Nayak
%
%  The finite element formulation of the equation to be solved,
%  using the method of weighted residuals and isoparametric element method
%  The natural triangular coordinate system is (0,0),(1,0) and (0,1)
%  Mapping functions are N1 = 1-U-V, N2 = U and N3 = V
%
%Problem: Parallel plate capacitor voltage distribution calculation.
%    Boundary conditions : Top plate voltage and Bottom plate voltage
%
%
%
%
clc; clear;
%%
%Coordinate to be analyzed
X_range = [0 1];
Y_range = [0 1];

dx = 0.2;
dy = 0.2;

nx = 1 + (X_range(2) - X_range(1)) / dx;
ny = 1 + (Y_range(2) - Y_range(1)) / dy;

%Meshgrid formation
[x,y] = meshgrid(X_range(1) : dx : X_range(2) , Y_range(1) : dy : Y_range(2));
D = delaunay(x,y);                             % Triangulation
triplot(D,x,y)

%Material parameter
disp('MATERIAL PARAMETER');
BC_top = input('Enter the voltage on Top plate = ');
BC_bottom = input('Enter the voltage on Bottom plate = ');
epsilon_0 = 8.85e-12;
epsilon_r = 5;
epsilon = epsilon_0 * epsilon_r;
Rho = 1e-9;

M = size(D,1);                           % Number of Elements
max_Node = max(max(D));                  % Total number of nodes
%Coordinate Separation
for i = 1 : M
    for j = 1 : 3
        xc{i}(j) = x( D(i , j) );        %Element wise coorinate storage
        yc{i}(j) = y( D(i , j) );
    end
end

syms U;
syms V;

%Element-wise processing
for i = 1 : M
    
    x1 = xc{i}(1);                          %Coordinate copy element-wise
    x2 = xc{i}(2);
    x3 = xc{i}(3);
    
    y1 = yc{i}(1);
    y2 = yc{i}(2);
    y3 = yc{i}(3);
    
    N1 = 1-U-V;                                    % Mapping functions
    N2 = U;
    N3 = V;
    N = [N1 N2 N3];                                % Mapping vector

    DelN = [ diff(N,U) ; diff(N,V) ];

    xv = [x1 ; x2 ; x3];
    yv = [y1 ; y2 ; y3];
    X = N * xv;
    Y = N * yv;

    %Jaccobian calculation
    
    J = [ diff( X , U ) diff( Y , U ) ; diff( X , V ) diff( Y , V ) ];
    
    k = inv(J) * DelN;
    k1 = epsilon .* k' * k * det(J);
    
    %Local stiffness matrix
    Ke{i} = vpaintegral(vpaintegral(k1 , V , [0 1-U]) , U , [0 1]);
    
    %Local forced vector calculation
    Fe{i} = vpaintegral(vpaintegral((Rho.* N' * det(J)) , V , [0 1-U]) , U , [0 1]);
    
end

%% Assembly of Global Stiffness Matrix
K = zeros(max_Node,max_Node); %Initializing global stiffness matrix
for p = 1:M                            %iterating over all elements
for i=1:3                              %iterating over all the rows of local stiffness matrix
for j= 1:3                             %iterating over all the columns of local stiffness matrix
K(D(p,i), D(p,j)) = K(D(p,i), D(p,j))+ Ke{p}(i,j);
%Assembling the values at global node by mapping values at local nodes
%(thru all elements) into it
end
end
end
%% Assembly of Global force matrix
Force= zeros(max_Node,1); %Initializing global force matrix
for p = 1:M
for i=1:3
Force(D(p,i)) = Force(D(p,i)) +Fe{p}(i);
end
end


%% Applying Boundary conditions
BCnodes_top = [];                       % intializing boundary nodes at top of parallel plate
BCnodes_bottom = [];                    % intializing boundary nodes at bottom of parallel plate
for i= 1:max_Node
if rem(i,ny)==0
BCnodes_top = [BCnodes_top, i];         % Boundary nodes at top electrode
end
if rem(i,ny)==1
BCnodes_bottom = [BCnodes_bottom, i];   % Boundary nodes at bottom electrode
end
end
for i = 1:length(BCnodes_top)           % loop over all top boundary nodes
nBC = BCnodes_top(i);                   % boundary node
K(nBC,:) = 0;                           % make the n-th row of global K matrix zero
K(nBC,nBC) = 1;                         % make the diagonal term of the global K matrix 1
Force(nBC) = BC_top;                    % assign the n-th row of Force to the known boundary value
end
for i = 1:length(BCnodes_bottom)        % loop over all bottom boundary nodes
nBC = BCnodes_bottom(i);                % boundary node
K(nBC,:) = 0;                           % make the n-th row of global K matrix zero
K(nBC,nBC) = 1;                         % make the diagonal term of the global K matrix 1
Force(nBC) = BC_bottom;                 % assign the n-th row of Force to the known boundary value
end
%% Solution to FEM
%Equation is of form [K][Phi]= [F]
Phi = (inv(K))*Force;                   %nodal voltage values
Phiout = ([reshape(Phi,[nx,ny] )]);     %rearanging nodal values in nx X ny matrix
%% Plotting solution
figure

pcolor(x,y,Phiout)
shading interp
colorbar
colormap 'jet'
title('Voltage distribution')
xlabel('distance in x')
ylabel('distance in y')