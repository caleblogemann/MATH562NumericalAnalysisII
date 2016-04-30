function [x0, k, r] = ConjugateGradient(A, b, tol, maxIter)
    x0 = zeros(size(b));
    k = 0;
    r0 = b - A*x0;
    r = norm(r0, inf);
    while(r(end) > tol && k < maxIter)
        k = k+1;
        Ar0 = A*r0;
        a = (r0'*Ar0)/(Ar0'*Ar0);
        x = x0 + a*r0;
        x0 = x;
        r0 = b - A*x0;
        r = [r, norm(r0,inf)];
    end
end
