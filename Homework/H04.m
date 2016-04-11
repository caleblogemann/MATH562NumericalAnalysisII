%% Problem 5
% part (a)
A = diag([-4, 2, 1, 1, 1]) + triu(rand(5,5),1);
v0 = ones(5, 1);
[v1, lam1, k1] = Pwr1(A, v0)
[v2, lam2, k2] = Pwr2(A, v0)

% part (b)
A = diag([9, 2, 1, 5, -8]) + triu(rand(5,5),1);
v0 = ones(5, 1);
[v1, lam1, k1] = Pwr1(A, v0)
[v2, lam2, k2] = Pwr2(A, v0)

%% Problem 6
A = diag([9, 2, 1, 5, -8]) + triu(rand(5,5),1);
v0 = ones(5, 1);
[v, lam, k] = Inv(A, v0, 8.8)

%% Problem 7
A = diag([9, 2, 1, 5, -8]) + triu(rand(5,5),1);
v0 = ones(5, 1);
[v, lam, k] = Ray(A, v0)

