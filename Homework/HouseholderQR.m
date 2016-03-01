function [W, R] = HouseholderQR(A)
    p = inputParser;
    p.addRequired('A', @Utils.isMatrix);
    p.parse(A);

    [m, n] = size(A);
    W = zeros(m, n);
    for k=1:n
        x = A(k:m, k);
        %vk = sign(x(1))*norm(x,2)*eye(m-k+1, 1) + x;
        vk = norm(x,2)*eye(m-k+1, 1) + x;
        vk = vk/norm(vk, 2);
        A(k:m, k:n) = A(k:m, k:n) - 2*vk*(vk'*A(k:m, k:n));
        W(k:m, k) = vk;
    end
    R = A;

