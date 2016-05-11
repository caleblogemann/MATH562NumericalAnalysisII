clear all
close all

format long

% fitting exp(sin(4t)) on [0 1]
% by polynomial of degree n-1
% with m points

m = 100; n = 15;
t = linspace(0,1,m)';

A=[]; %Vandermonde matrix
for i=1:n
    A = [A t.^(i-1)];
end

b = exp(sin(4*t));
%xx = A\b; 
%b = b/xx(n);
b=b/2006.787453080206; % normalize so that x(n) = 1

%%
% Solve the least square problem; assume to be accurate
% to obtain the parameters: Kappa(A), Theta, Eta. 
['Computing Kappa, Theta, Eta']
x = A\b; y = A*x;
kappa = cond(A);
theta = asin(norm(b-y)/norm(b));
eta = norm(A)*norm(x)/norm(y);
kappa
theta 
eta


%% Householder Triangularization without column pivoting
% explicit Q 
['Householder, No Column Pivoting, Explicit Q']
[Q R] = qr(A);
x = R\(Q'*b);
x(n)
norm(x(n)-1)
norm(x(n)-1)/kappa



%% Householder Triangularization without column pivoting of Augmented system
% implicit Q
['Householder, No Column Pivoting, Implicit Q']
[Q R] = qr([A b]);
Qb = R(1:n,n+1);
R = R(1:n,1:n);
x = R\Qb;
x(n)
norm(x(n)-1)
norm(x(n)-1)/kappa


%% Householder QR with column pivoting
['Householder, Column Pivoting']
x = A\b;
x(n)
norm(x(n)-1)
norm(x(n)-1)/kappa

%% Modified Gram-Schmidt Orthogonalization
% explicit Q
['Modified Gram-Schmidt, Explicit Q']
[Q R]=mgs(A);
x = R\(Q'*b);
x(n)
norm(x(n)-1)
norm(x(n)-1)/kappa

%an example that mgs loses orthogonality
AA = [0.7 0.70711; 0.70001 0.70711];
[Q R] = qr(AA); 
norm(Q'*Q-eye(2))
[Q R] = mgs(AA);
norm(Q'*Q-eye(2))


%% Modified Gram-Schmidt Orthogonalization of Augmented system
% implicit Q
['Modified Gram-Schmidt, Imlicit Q']
[Q R] = mgs([A b]);
Qb = R(1:n,n+1);
R = R(1:n,1:n);
x = R\Qb;
x(n)
norm(x(n)-1)
norm(x(n)-1)/kappa



%% Normal Equation 
['Normal Equations']
x = (A'*A)\(A'*b);
x(n)
norm(x(n)-1)
norm(x(n)-1)/(kappa)
norm(x(n)-1)/(kappa^2)



%% SVD
['SVD']
[U S V] = svd(A);
x = V*(S\(U'*b));
x(n)
norm(x(n)-1)
norm(x(n)-1)/kappa

