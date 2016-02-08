function [U, S, V] = SVD2D(A)
    % create possible unit vectors v1
    theta = 0:.01:2*pi;
    v = [cos(theta); sin(theta)];

    Av = A*v;
    normAv = arrayfun(@(i) norm(Av(:, i), 2), 1:length(theta));
    [s1, i] = max(normAv);

    % v_i maximizes 2-norm of Av
    v1 = v(:, i);
    % s1*u1 = A*v1
    Av1 = A*v1;
    s1 = norm(Av1,2);
    u1 = Av1/s1;
    
    % v2 should be orthogonal to v1
    % generate random vector
    w = rand(2,1);
    v2 = w - (w'*v1)*v1;
    v2 = v2/norm(v2);

    % s2*u2 is A*v2
    Av2 = A*v2;
    s2 = norm(Av2,2);
    u2 = Av2/s2;

    U = [u1, u2];
    V = [v1, v2];
    S = diag([s1, s2]);

    % now plot 
    circle_image(A);
    % plot v vectors
    arrow([0, 0], v1, 'b');
    arrow([0, 0], v2, 'b');
    % plot u vectors
    arrow([0, 0], s1*u1, 'r');
    arrow([0, 0], s2*u2, 'r');
end
