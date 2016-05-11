function [lambda q] = RayleighQuotientIteration(A,v0,maxIter)

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

%q = randn(m,1); q=q/norm(q);
q=v0; q = q/norm(q);
lambda = q'*A*q;

for k=1:maxIter
    %lambda0 = lambda;
    if cond(A-lambda*eye(m))>1e12
        break
    end
    w = (A-lambda*eye(m))\q;
    q = w/norm(w);
    lambda = q'*A*q;
    %if abs(lambda0-lambda)<1e-9
    %    break
    %end
end
