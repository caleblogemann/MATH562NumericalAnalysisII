function [lambda q] = PowerIteration(A,maxIter)

if ~isreal(A)
    ['Input is not real matrix\n']
    return
end

%if ~issymmetric(A)
%    ['Input is not symmteric matrix\n']
%    return
%end


[m n]=size(A);
if m~=n 
    ['Input is not square matrix']
    return
end

q = randn(m,1); q=q/norm(q);
lambda = 0;

for k=1:maxIter
    w = A*q;
    q = w/norm(w);
    lambda = q'*A*q;
end