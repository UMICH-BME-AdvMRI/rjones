function Rx=xrot(phi)
% function Rx=xrot(phi)
%
% Generate rotation matrix for rotation phi about x-axis
%
% INPUTS: phi = rotation angle (radians)
% OUTPUTS: Rx = 3x3 rotation matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% PLEASE NOTE
%   This function was from Brian Hargreaves Bloch Equation Simulation webpage
%   http://mrsrl.stanford.edu/~brian/bloch/



Rx = [1 0 0; 0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

