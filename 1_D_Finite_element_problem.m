%This is a MATLAB programme to solve the 1-D problem
%NAME: Sandeep Mohan Nayak
%
%
%  finite element method to solve the following equation.
%  The finite element formulation of the equation to be solved,
%  using the method of weighted residuals
%
%        d^2u
%    *********** - u = -x, 0 < x < 1
%        dx^2
%
% BOUNDARY CONDITION u(0) = 0 , u(l) = 0
%

clc;
clear;
close all;

%Taking the input data
%Pre-processing
len = 'Enter the length\n';
l = input(len);
node = 'Enter the number of nodes\n';
n = input(node);

e = n-1;                              %Number of Elements
dx = l/e;                             %Step size
A = 0:dx:l;                           %Coordinate value
Al = length(A);


%**************************************************************************
%Matrix to take the coordinate value of the elements
%CM = [Material_type.......Element_type........Begin_Node.........End_Node]
%


CM = zeros(e,4);                     % Element connectivity data matrix
Kg = zeros(n);                       % Global stiffness matrix
Fg = zeros(n,1);                     % Global forced vector


%Element connectivity data matrix

for i = 1:Al-1
    CM(i) = 1;                       % Material Type
    CM(i+e) = 1;                     % Element type
    CM(i+(2*e)) = A(i);              % Begin node
    CM(i+(3*e)) = A(i+1);            % End node
end

%Boundary condition information
nd1 = 1;                             % Node number - 1
nd2 = n;                             % Node number - 2
Bc(1) = 0;                           % Boundary Condition @ Node-1
Bc(2) = 0;                           % Boundary Condition @ Node-2


%%
%**************************************************************************
% Processing...


a = 0;
b = 0;
c = 0;
for j = 1:Al-1
    
syms x;
%Capturing the actual node coordinate
xc = CM(j+(2*e));                    % Starting node position
xp = CM(j+(3*e));                    % Ending node position


H1 = (xp - x)/(xp - xc);
H2 = (x - xc)/(xp - xc);
H = [H1 ; H2];
K = (-diff(H)*diff(H'))-(H*H');


Ke = int(K,x,[xc xp]);               % Local stiffness matrix
Fe = int(-H*x,x,[xc xp]);            % Local forced vector


%*************************************************************************


%Global stiffness matrix formation Kg
Kg(j+a) = Ke(1) + b;                 %Diagonal Number
Kg(j+1+a) = Ke(2);
Kg(j+n+a) = Ke(3);                   
Kg(j+n+1+a) = Ke(4);                 %Diagonal Number
a = a + n;
b = Ke(4);


%*************************************************************************
%Global forced vector Formation Fg


Fg(j) = Fe(1) + c;
Fg(j+1) = Fe(2);
c = Fe(2);
end


%Implementation of Algorithm-1 to Modify the Global stiffness matrix


Kg1 = Kg;
Kg1(nd1,:) = 0;
Kg1(nd2,:) = 0;
Kg1(nd1,nd1) = 1;
Kg1(nd2,nd2) = 1;


%Implementation of Algorithm-1 to Modify the Global Forced matrix


Fg1 = Fg;                            % Modified Stiffness matrix ac to Bcs
Fg1(nd1) = Bc(1);
Fg1(nd2) = Bc(2);

%Solutiion
U1 = inv(Kg1)*Fg1;


%**************************************************************************
%Post Processing........


plot(A,U1)
xlabel('l')
ylabel('u')
title('Finite Element Solution of 1-D problem')