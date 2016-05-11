clear all
close all
%format long
format short

%given matrix (symmteric)
m = 6; 
%A = randn(m); 
%A = invhilb(m);
A = full(gallery('tridiag',m,-1,2,-1));
%A = magic(m);
A = A'*A+rand(m); 
A0=A;

% if issymmetric(A) ~= 1
%     ['A is not symmteric']
%     return
% end

[QE Lam] = eig(A);   % "exact" e-value
%Lam = diag(Lam);
Lam = sort(diag(Lam),'descend');
QE = QE(:,[m:-1:1]);


maxIter = 100;

%%%%%%%%%
%PURE QR
%%%%%%%%%
%phase I: householder to Hessenberg
A = hess(A);

%phase II: "Pure" QR
for k=1:maxIter
    [Q R] = qr(A);
    A = R*Q;
end

lam = diag(A);

%[Lam lam]


%return

%%%%%%%%%%%%%%
%Simultaneous
%Iteration
%%%%%%%%%%%%%%
A=A0;

%Phase I:
A = hess(A);

%Phase II:
%[Q2 R0] = qr(randn(m));
Q2=eye(m);
for k=1:maxIter
    Z = A*Q2;
    [Q2 R2] = qr(Z);
end

lam2 = diag(R2);

%return 

%%%%%%%%%%%%%%
%Simultaneous
%inverse
%Iteration
%%%%%%%%%%%%%%
A=A0;

%Phase I:
A = hess(A);

%Phase II:
%[Q2 R0] = qr(randn(m));
Q3=flip(eye(m));
for k=1:maxIter
    Z = A\Q3;
    [Q3 R3] = qr(Z);
end

lam3 = flip(1./diag(R3));



%%%%%%%%%%%%%%
%QR
%with shifts
%Rayleigh quotien shift
%%%%%%%%%%%%%%
A=A0;

%Phase I:
A = hess(A);

%Phase II
for k=1:maxIter
    mu = A(m,m);
    [Q4 R4] = qr(A-mu*eye(m));
    A = R4*Q4+mu*eye(m); 
    
    %A 
    %pause
end

%lam4 = sort(diag(A),'descend');
lam4 = diag(A);

%return 
%%%%%%%%%%%%%%
%Practical QR
%with shifts
%Rayleigh quotien shift
%%%%%%%%%%%%%%
['Practical QR']
A=A0;

%Phase I:
A = hess(A);

U = eye(m);
%k = 0; 
for n = m:-1:2
    k = 0;
    %while abs(A(n,n-1))>1e-30 & k<maxIter
    mu = 100; oldmu = 1;
    while abs(oldmu-mu)>1e-10
        oldmu = mu;
        k=k+1;
        mu = A(n,n);
        [Q5 R5] = qr(A-mu*eye(m));
        A = R5*Q5+mu*eye(m);
        U = U*Q5;
    end
    A(n,n-1) = 0; A(n-1,n) = 0;
        
    A 
    pause
end
 
lam5 = sort(diag(A),'descend');



[Lam lam lam2 lam3 lam4 lam5]
%[abs(lam-Lam) abs(lam2-Lam) abs(lam3-Lam) abs(lam4-Lam)]