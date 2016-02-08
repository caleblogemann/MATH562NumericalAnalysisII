function circle_image(A)

% Script File : circle_image.m
%
% This script sets up a (discrete) parameterization of the unit
% circle and computes the image of it under the map define by
% matrix multiplication by a 2x2 matrix A.  It assumes A has
% been previous defined.  Both the circle and its image are
% plotted in the same figure with image points of (1,0) and
% (0,1) indicated by *'s.  

t = linspace(0,2*pi,361);
c = cos(t);
s = sin(t);
plot(c,s, 'b', 'Linewidth', 2)
axis equal
hold on
x = [c;s];
y = A*x;
plot(y(1,:),y(2,:),'r', 'Linewidth', 2)
