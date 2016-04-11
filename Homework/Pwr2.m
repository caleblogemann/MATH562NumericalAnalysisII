function [v, lam, k] = Pwr2(A, v0)

    % normalize v0
    v = v0/norm(v0,2);
    % initiate lam variables so while loop condition works for first iteration
    lam0 = 0;
    lam = 1;
    k = 0;
    while k < 500 && abs(lam - lam0) > 1e-8
        % update next step v0 acts like v^(k-1)
        k = k + 1;
        v0 = v;
        lam0 = lam;

        w = A*v0;
        % normalize
        v = w/norm(w,2);
        % Rayleigh quotient
        lam = v'*A*v;
    end
end
