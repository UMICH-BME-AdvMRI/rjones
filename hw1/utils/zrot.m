function Rz=zrot(phi)
% function Rz=zrot(phi)
%
% Generate rotation matrix for rotation phi about z-axis
%
% INPUTS: phi = rotation angle (radians)
% OUTPUTS: Rz = 3x3 rotation matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% PLEASE NOTE
%   This function was from Brian Hargreaves Bloch Equation Simulation webpage
%   http://mrsrl.stanford.edu/~brian/bloch/

Rz = [cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0; 0 0 1];

