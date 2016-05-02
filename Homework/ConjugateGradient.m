function [x0, k, r] = ConjugateGradient(A, b, tol, maxIter)
    x0 = zeros(size(b));
    k = 0;
    r0 = b;
    p0 = r0;
    r = norm(r0, inf);
    while(r(end) > tol && k < maxIter)
        a = (r0'*r0)/(r0'*A*p0);
        x = x0 + a*p0;
        r1 = r0 - a*A*p0;
        bn = (r1'*r1)/(r0'*r0);
        p0 = r1 + bn*p0;

        % move to next step
        k = k+1;
        x0 = x;
        r0 = r1;
        r = [r, norm(r0,inf)];
    end
end
