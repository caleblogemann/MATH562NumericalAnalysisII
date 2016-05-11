clear all
close all

format long

ep = 1e-10;
A = [1 1 1; ep 0 0; 0 ep 0; 0 0 ep];

[Q1 R1] = mgs(A);
[Q2 R2] = gs(A);

Q1

Q2


Q1(:,2).'*Q1(:,3)

Q2(:,2).'*Q2(:,3)