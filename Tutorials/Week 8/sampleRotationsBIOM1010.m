function sampleRotationsBIOM1010(u,theta)

% Sample rotations of reference frame axes about axis u by angle theta
% measured in (rad). u does not need to be normalised.

% Rotate from xy plane to yz plane
u=u/norm(u); % This is the bit missing in the lecture which is why it didn't work!

q = [cos(theta/2) sin(theta/2)*u]

x = [1 0 0]
xdash = quatrot(x,q)'

y = [0 1 0]
ydash = quatrot(y,q)'

z = [0 0 1]
zdash = quatrot(z,q)'

end

function vrot = quatrot(v,q)

% Use quaternion to rotate vector
% 
% v: vector to rotate
% q: unit quaternion to perform rotation around vector (x,y,z) using 
% the quaternion q = w + xi +yj +zk
%
% Ensure column vectors
q = q(:);
v = v(:);

% Check size
if length(q)~=4
    error('Quaternion must contain four elements.');
end

% Ensure q is unit quaternion
normq = sqrt(sum(q.^2));
q = q./normq;

% vrot = qvq*
vrot = quatmultiply(quatmultiply(q,[0; v]),quatinv(q));
vrot = vrot(2:4);

end

function s = quatmultiply(p,q)

% Multiply two quaternions
%
% Input
% p: quaternion p(1) + p(2)i + p(3)j + p(4)k
% q: quaternion q(1) + q(2)i + q(3)j + q(4)k
%
% Output
% s: quaternion product of p and q

% Ensure column vector
q = q(:);
p = p(:);

% Check size
if length(q)~=4 || length(p)~=4
    error('Quaternion must contain four elements.');
end

s =  [p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4)...
    ; p(1)*q(2) + p(2)*q(1) + p(3)*q(4) - p(4)*q(3)...
    ; p(1)*q(3) - p(2)*q(4) + p(3)*q(1) + p(4)*q(2)...
    ; p(1)*q(4) + p(2)*q(3) - p(3)*q(2) + p(4)*q(1)];

end

function qinv = quatinv(q)

% Calculates inverse of quaternion by negating ijk components
%
% Input::
% v    : vector it needs to be rotated around
% theta: angle to be rotated by
%
% Output::
% q : quaternion output, with vector component normalised and angle halved


% Ensure column vector
q = q(:);

% Check size
if length(q)~=4
    error('Quaternion must contain four elements.');
end

qinv = [q(1); -q(2:4)];

end
