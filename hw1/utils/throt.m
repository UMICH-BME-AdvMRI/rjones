function Rth=throt(phi,theta)
% function Rth=throt(phi,theta)
%
% Generate rotation matrix for rotation about arbitrary transverse (xy)
% axis; rotation of phi radians about the axis defined by y=x*tan(theta)
%
% INPUTS: phi = rotation angle (radians), 
%         theta = to define rotational axis y=x*tan(theta) (radians) 
% OUTPUTS: Rth = 3x3 rotation matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% PLEASE NOTE
%   This function was from Brian Hargreaves Bloch Equation Simulation webpage
%   http://mrsrl.stanford.edu/~brian/bloch/

Rz = zrot(-theta);
Rx = xrot(phi);
Rth = inv(Rz)*Rx*Rz;

