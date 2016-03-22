inverseHilbert = invhilb(6);
A = inverseHilbert(:, 1:5);
[m, n] = size(A);
b = [463; -13860; 97020; -258720; 291060; -116424];
xExact = [1; 1/2; 1/3; 1/4; 1/5];
y = A*x;
theta = acos(norm(y)/norm(b));
eta = norm(A)*norm(xExact)/norm(y);
kappa = cond(A); 

condby = 1/cos(theta)
condbx = kappa/(eta*cos(theta))
condAy = kappa/cos(theta)
condAx = kappa + kappa^2*tan(theta)/eta

% Householder QR
[Q, R] = qr(A, 0);
x = R\(Q'*b);
E = norm(x - xExact);
M = {'Householder QR'};

% Householder QR of augmented matrix
[Q, R] = qr([A, b], 0);
Qb = R(1:n, n+1);
R = R(1:n, 1:n);
x = R\Qb;
E = [E; norm(x - xExact)];
M = [M, {'Householder QR of augmented matrix'}];

% Modified Gram-Schmidt QR
[Q, R] = mgs(A);
x = R\(Q'*b);
E = [E; norm(x - xExact)];
M = [M, {'Modified Gram-Schmidt QR'}];

% Modified Gram-Schmidt QR of augmented matrix
[Q, R] = mgs([A, b]);
Qb = R(1:n, n+1);
R = R(1:n, 1:n);
x = R\Qb;
E = [E; norm(x - xExact)];
M = [M, {'Modified Gram-Schmidt QR of augmented matrix'}];

% Normal Equation
x = (A'*A)\(A'*b);
E = [E; norm(x - xExact)];
M = [M, {'Normal Equations'}];

% SVD
[U, S, V] = svd(A, 0);
x = V*(S\(U'*b));
E = [E; norm(x - xExact)];
M = [M, {'SVD'}];

table(E, 'VariableNames', {'Error'}, 'RowNames', M)
