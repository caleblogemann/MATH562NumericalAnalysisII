function [Q R] = mgs(A)

if nargin ~= 1,
    error('You need to imput only one argument.');
end

[m n]=size(A);

for j= 1:n
    R(j,j)=norm(A(:,j));
    A(:,j)=A(:,j)/R(j,j);
    R(j,j+1:n)=A(:,j)'*A(:,j+1:n);
    A(:,j+1:n)=A(:,j+1:n)-A(:,j)*R(j,j+1:n);
end

Q = A;

return