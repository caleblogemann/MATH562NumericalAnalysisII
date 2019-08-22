A = [1, 2, 3; 4, 5, 6; 7, 8, 7; 4, 2, 3; 4, 2, 2];
[W, R] = HouseholderQR(A)
[m, n] = size(A);
q1 = eye(m,1);
q2 = zeros(m,1);
q2(2) = 1;
q3 = zeros(m,1);
q3(3) = 1;
q4 = zeros(m,1);
q4(4) = 1;
for k = n:-1:1
    vk = W(k:m, k);
    q1(k:m) = q1(k:m) - 2*vk*(vk'*q1(k:m));
    q2(k:m) = q2(k:m) - 2*vk*(vk'*q2(k:m));
    q3(k:m) = q3(k:m) - 2*vk*(vk'*q3(k:m));
    q4(k:m) = q4(k:m) - 2*vk*(vk'*q4(k:m));
end
Q = [q1, q2, q3, q4];
Q*Q'

A = [.5, -2.5; -.5, 2.5; -.5, -7.5; -.5, -7.5];
[W, R] = HouseholderQR(A)

[m, n] = size(A);
q1 = eye(m,1);
q2 = zeros(m,1);
q2(2) = 1;
q3 = zeros(m,1);
q3(3) = 1;
q4 = zeros(m,1);
q4(4) = 1;
for k = n:-1:1
    vk = W(k:m, k);
    q1(k:m) = q1(k:m) - 2*vk*(vk'*q1(k:m));
    q2(k:m) = q2(k:m) - 2*vk*(vk'*q2(k:m));
    q3(k:m) = q3(k:m) - 2*vk*(vk'*q3(k:m));
    q4(k:m) = q4(k:m) - 2*vk*(vk'*q4(k:m));
end
Q = [q1, q2, q3, q4];
