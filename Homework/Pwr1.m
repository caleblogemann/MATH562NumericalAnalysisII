function [v, lam, k] = Pwr1(A, v0)

    % normalize v0 
    v = v0/norm(v0,2);
    % replace v0 so that while loop condition works
    v0 = zeros(size(v));
    % initiate v so while loop doesn't quite at first check
    k = 0;
    lam = 0;
    while k < 500 && norm(v - v0, 2) > 1e-8
        % update next step v0 acts like v^(k-1)
        k = k + 1;
        v0 = v;

        w = A*v0;
        % normalize
        v = w/norm(w,2);
        % Rayleigh quotient
        lam = v'*A*v;
    end
end
