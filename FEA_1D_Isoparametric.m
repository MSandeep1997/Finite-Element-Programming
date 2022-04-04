%This is a MATLAB programme to solve the 1-D FEM programme using
%Isoparametric elements
%NAME: Sandeep Mohan Nayak
%
%
%  finite element method to solve the following equation.
%  The finite element formulation of the equation to be solved,
%  using the method of weighted residuals and isoparametric element method
%  The natural coordinate system is -1 to 1
%  Mapping functions are N1 = (1-u)/2 and N2 = (1+u)/2
%
%        d^u
%    *********** - u = -x, 0 < x < 1
%        dx^2
%
% BOUNDARY CONDITION u(0) = 0 , u(l) = 0
%
clc;
clear;
close all;

% Taking input data
len = 'enter the length\n';
l = input(len);
nod = 'enter the number of nodes\n';
n = input(nod);

dx = l/n;
cod = 0:dx:l;                     % Coordinates in the vector
lc = length(cod);

syms u;
%Processing
for i=1:lc-1
    x1 = cod(i);                   % Coordinates assign to a variable
    x2 = cod(i+1);                 % Coordinates assign to another variable
    xv = [x1 ; x2];
    
    N1 = (1-u)/2;                  % Mapping functions
    N2 = (1+u)/2;
    
    N = [N1 N2];
    Nd = diff(N);
    
    x = N * xv;
    J = diff(x);                    % Jaccobian Calculation
    
    
    k = Nd * inv(J);
    
    K = (k' * k + N' * N) * J;
    F = N' * N * xv * J;

    Kc{i} = int(K,u,[-1 1]);        % Local stiffness matrix
    Fc{i} = int(F,u,[-1,1]);        % Local forced vector
end

%Formation of global stiffness matrix
for j=1:lc-1
    if j == 1
        a = [0 0;0 0];
        Kn(j:j+1,j:j+1) = Kc{j}+a;
    else
        a = [Kc{j-1}(2,2) 0;0 0];
        Kn(j:j+1,j:j+1) = Kc{j}+a;
    end
end
%Formation of global forced vector
for j=1:lc-1
    if j == 1
        b = [0;0];
        Fn(j:j+1) = Fc{j}+b;
    else
        b = [Fc{j-1}(2);0];
        Fn(j:j+1) = Fc{j}+b;
    end
end
Fn = Fn';

%Assignment of Boundary conditions
Kn1 = Kn;
Kn1(1,:) = 0;
Kn1(lc,:) = 0;
Kn1(1,1) = 1;
Kn1(lc,lc) = 1;

Fn1 = Fn;
Fn1(1) = 0;
Fn1(lc) = 0;

%%Post processing

U = inv(Kn1)*Fn1;

plot(cod,U)
xlabel('l')
ylabel('u')
title('Finite Element Solution of 1-D problem')