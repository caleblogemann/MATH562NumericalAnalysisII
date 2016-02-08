function arrow(start,finish,s1)

% Draws an arrow from a starting point to a finish point
% in the current (2d) plot.  The arrow head will be given the
% color indicated by the character string s1 ('r','b', etc.).

% make sure both vectors have the same orientation.
% and check to see if finish - start defines a non-zero vector.

if (size(start) ~= size(finish))
  finish = finish';
end
if (norm(finish - start) == 0)  % do nothing if the vector = 0
  return
end

% the angle between the two sides of arrow will
% 30 degrees = 2*15 degrees

theta=pi/12; 
cos_theta=cos(theta); 
sin_theta=sin(theta);

x1 = start(1); y1 = start(2);      % extract data
x2 = finish(1); y2 = finish(2);

plot([x1,x2],[y1,y2],'k');         % plot line seqment

% now make the arrow head

v = axis;   % get scaling for current plot

% size of arrow head will be 1/50 of the diagonal 
% of the whole picture 

len = sqrt((v(2) - v(1))^2 + (v(4) - v(3))^2)/50;

p = zeros(4,1); q = zeros(4,1);  % vertices of arrow head
p(1) = x2; q(1) = y2;
p(4) = x2; q(4) = y2;

% define a unit vector from finish to start
% a rotations through +theta,-theta radians

u = [x1 - x2;y1 - y2]/norm(finish-start); 
Q = [cos_theta, -sin_theta; sin_theta, cos_theta];

edge_vec = len*Q*u;
p(2) = p(1) + edge_vec(1);
q(2) = q(1) + edge_vec(2);

Q = [cos_theta, sin_theta; -sin_theta, cos_theta];

edge_vec = len*Q*u;
p(3) = p(1) + edge_vec(1);
q(3) = q(1) + edge_vec(2);

fill(p,q,s1)
