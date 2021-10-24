function [AHRS] = mainAHRS(Acc,Mag,Gyro,timeInput)

% Attitude and heading reference system (AHRS) for estimating device 
% orienation and acceleration components resolved in the global frame. 
% Requires triaxial acceleration, triaxial magnetometer, and triaxial
% gryroscope.
%
% Input::
% Acc:       Nx3 matrix of acceleration components measured by the device 
%            in units of m/s^2.
% Mag:       Nx3 matrix of magnetic components measured by the device; 
%            units not too important but coming close to a ferromagnetic 
%            metal can distort earth's magnetic field. Calibration for soft- 
%            and hard-iron effects is very important.
% Gyro:      Nx3 matrix of gyroscope components measured by the device 
%            in units of rad/s. Calibration is important.
% timeInput: Either a scaler giving sampling frequency, or an array of 
%            sample times in seconds.       
%
% Output::
% AHRS.q     Nx3 quaternion orientations at each frame.
% AHRS.R:    3x3xN rotation matrix defining orienation of device in global
%            frame.
% AHRS.Acc:  Nx3 matrix of acceleration x,y,z acceleration components
%            resolved in the global frame.
% AHRS.Mag:  Nx3 matrix of magnetometer values, resolved in global frame. 

% Author: Stephen Redmond, 10 Oct 2014. Version 1.

% Number of time samples
N = size(Acc,1);

% Sampling times
if ~all(size(timeInput)== [1 1])
    t = timeInput; % Possibly non-uniform sampling. Sample times given.
else
    t = ones(N,1)/timeInput; % Uniform sampling
end

% Initialise memory
qGlobal(N,4)=NaN;
RGlobal(1:3,1:3,1:N) = NaN;
AccGlobal(1:N,1:3) = NaN;
MagGlobal(1:N,1:3) = NaN;

% Initialise orientations
qGlobal(1,:) = [1 0 0 0];
RGlobal(1:3,1:3,1) = eye(3);

% Initialise global accelerations using best guess of orientation
AccGlobal(1,1:3) = (RGlobal(1:3,1:3,1)*Acc(1,1:3)')';

initialisationDuration = 1; % Time spent initialising in seconds

for i = 2:N

    dt = t(i) - t(i-1);
    
    % Gear shifting
    if t(i)>initialisationDuration % After initialisation is over
        % muXXX is the fractional rotation tuning parameter.
        % Defines orientation responsiveness to gravity or magnetometer changes. 
        % muAccGyro=0 ignores acc and uses on gyro. 
        % muAccGyro=1 ignores gyro and uses only acc.
        % muMag = 0 ignores magnetometer and 1 trusts it completely 
        
        %%% User defined %%%
        muMag = 0.01*dt; % Correction rate for magnetometer alignment.
        muAcc = 0.2*dt; % Correction rate for accelerometer alignment.

        %%% Removed this gear-shifting for BIOM1010 class
%         if abs(norm(Acc(i,:))-9.81)<1
%             %%% User defined %%%
%             muAcc = 4.0*dt; % Correction rate for accelerometer alignment 
%         else
%             %%% User defined %%%
%             muAcc = 0.2*dt; % Correction rate for accelerometer alignment 
%         end

    else
        % Fast initialisation at start
        muAcc = 8*dt; 
        muMag = 8*dt;        
    end
        
    %%% Correct for gyros
        
    % Get quaternion from gyro in device frame and rotate
    qAngularVelocityDevice = angularVelocityToQuaternion(Gyro(i,:),dt);
    % Convert to back to global frame
    qGlobal(i,:) = quatmultiply(qGlobal(i-1,:),qAngularVelocityDevice);
%     qGlobal(i,:) = qGlobal(i-1,:);
    
    %%% Correct for accelerometer
    
    % Get next rotation required to align accelerometer/gravity with 'up'.
    % muAcc is the fraction of the way we move towards the gravity target. 
    % muAcc=1 means we trust accelerometer.
    
    % Convert acceleration to global frame
    AccGlobalTemp = quatrot(Acc(i,:),qGlobal(i,:));
    % Get fraction of rotation from acc in global to up
    [qRotTowardsUp] = getRotationQuaternion(AccGlobalTemp(:)',[0 0 1],muAcc);
    % Rotate global frame towards 'up'
    qGlobal(i,:) = quatmultiply(qRotTowardsUp,qGlobal(i,:));
    
    %%% Correct for magnetometer
    
    % Get next rotation around vertical to align measured global xy 
    % components of magnetic measurment with north.
    % muMag is the fraction of the way we move towards the mag target. 
    % muMag=1 means we trust magnetometer.
    
    % Transform magnetometer reading into global frame
    MagGlobalTemp = quatrot(Mag(i,:),qGlobal(i,:));
    % Get fraction of rotation from mag in global to north (xy components only)
    [qRotTowardsNorth] = getRotationQuaternion([MagGlobalTemp(1:2)' 0],[1 0 0],muMag);
    % Rotate global frame towards 'up'
    qGlobal(i,:) = quatmultiply(qRotTowardsNorth,qGlobal(i,:));
    
    %%%
   
    RGlobal(1:3,1:3,i) = getRotationMatrixFromQuat(qGlobal(i,:)); % Convert to rotation matrix

    % Project acceleration components into global frame one last time
    AccGlobal(i,1:3) = (RGlobal(1:3,1:3,i)*Acc(i,1:3)')';
    MagGlobal(i,1:3) = (RGlobal(1:3,1:3,i)*Mag(i,1:3)')';
    
end

AHRS.t = t;
AHRS.q = qGlobal;
AHRS.R = RGlobal;
AHRS.Acc = AccGlobal;
AHRS.Mag = MagGlobal;

end

function [qRot] = getRotationQuaternion(sensorVector,Ref,mu)
% Calculate quaternion rotation between where the device accleration vector
% is pointing (where the device believes down is) and down in the global
% frame
% 
% Input:: 
% sensorVector: 1x3 array of xyz accelerations in the device frame
%
% Output::  
% qAcc:         Single quaternion defining rotation required to get device
%               acceleration vector aligned with the global up

if mu<0
    mu = 0;
elseif mu>1
    mu = 1;
end

v1=sensorVector;
v2=Ref; % Reference vector in the global frame 
        % For gravity use [0 0 1](points up along z!)
        % For magnetic north use [1 0 0]
axisVector =cross(v1,v2);
axisNorm = sqrt(sum(axisVector.^2,2));
if axisNorm ==0
    qRot = [1 0 0 0];
else
    axisVectorNorm =  bsxfun(@rdivide,axisVector,axisNorm);
    % find the angle between the two vectors
    dotProductVectors=dot(v1',v2');
    magV1=sqrt(sum(v1.^2,2));
    magV2=sqrt(sum(v2.^2,2));
    angle=mu*acos(dotProductVectors'./(magV1.*magV2));
    q0 = cos(angle./2);                          % q.w
    q1 = axisVectorNorm(:,1).*sin(angle./2);     % q.x
    q2 = axisVectorNorm(:,2).*sin(angle./2);     % q.y
    q3 = axisVectorNorm(:,3).*sin(angle./2);     % q.z
    q = [q0 q1 q2 q3];
    % normalise to unit quaternions
    qNorm = sqrt(sum(q.^2,2));
    qRot = bsxfun(@rdivide,q,qNorm);
end

end

function q = angularVelocityToQuaternion(w,dt)

% Convert angular velocity vector into a quaternion for rotation. dt is the
% time that the rotation is performed for (assuming the angular velocity is 
% constant over that period) - usually this will be the sampling interval.
%
% Input::
% w:    angular velocity vector (angular velocities around three orthogonal 
% axes) in units of radians per second.
% dt:   time interval to rotate for, which determines the total angle of
% rotation
%
% Output::
% q:    quaternion which performs rotation through required angle around
% axis defined by angular velocity vector

angularSpeed = norm(w);  % Modulus of angular velocity (angular speed)
angle = angularSpeed*dt; % Total angle to rotate by
w = w/angularSpeed; %Normalise angular velocity vector for quaternion

q(1) = cos(angle/2);
q(2) = w(1)*sin(angle/2);
q(3) = w(2)*sin(angle/2);
q(4) = w(3)*sin(angle/2);

end

function A = getRotationMatrixFromQuat(q)

% Convert quaternion to rotation matrix
% From: http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
% 
% Input:: 
% q: Quaternion corresponding to rotation matrix A
% 
% Input::
% A: Rotation matrix derived from quaternion q

w = q(1); x = q(2); y = q(3); z = q(4);

n = w*w + x*x + y*y + z*z;

if n == 0
    s = 0; 
else
    s = 2 / n;
end

wx = s*w*x; 
wy = s*w*y; 
wz = s*w*z;

xx = s*x*x; 
xy = s*x*y; 
xz = s*x*z;

yy = s*y*y; 
yz = s*y*z; 
zz = s*z*z;

A = [ 1 - (yy + zz),         xy - wz,          xz + wy  ;...
           xy + wz ,    1 - (xx + zz),         yz - wx  ;...
           xz - wy ,         yz + wx,     1 - (xx + yy) ];
       
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


