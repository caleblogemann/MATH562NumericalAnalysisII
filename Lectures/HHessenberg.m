function [A dv] = HHessenberg(A)

[m n] = size(A);
if m~=n 
    ['Input is not a square matrix']
    return
end

dv = zeros(m-2,1);

for k=1:m-2
    x = A(k+1:m,k);
    v = x; 
    if x(1) == 0 
        s = 1;
    else
        s = sign(x(1));
    end
    v(1) = v(1) + s*norm(x);
    v = v/norm(v);
    A(k+1:m,k:m) = A(k+1:m,k:m) - 2*v*(v'*A(k+1:m,k:m));
    A(1:m,k+1:m) = A(1:m,k+1:m) - 2*(A(1:m,k+1:m)*v)*v';
    
    %if want to save v in the subdiagonal of A, uncomment these two lines.
    %A(k+2:m,k) = v(2:end);
    %dv(k) = v(1);
end

