function [v, lam, k] = Ray(A, v0)
    % set up
    v = v0/norm(v0,2);;
    lam = v'*A*v;
    % make sure while condition works initially
    lam0 = 0;
    k = 0;
    m = size(A);
    while k < 500 && abs(lam - lam0) > 1e-8
        % update next step v0 acts like v^(k-1)
        k = k + 1;
        v0 = v;
        lam0 = lam;

        % inverse iteration
        w = (A - lam0*eye(m))\v0;
        % normalize
        v = w/norm(w,2);
        % Rayleigh quotient
        lam = v'*A*v;
    end
end
