function Ry=yrot(phi)
% function Ry=yrot(phi)
%
% Generate rotation matrix for rotation phi about y-axis
%
% INPUTS: phi = rotation angle (radians)
% OUTPUTS: Ry = 3x3 rotation matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Brian Hargreaves 
%
% PLEASE NOTE
%   This function was from Brian Hargreaves Bloch Equation Simulation webpage
%   http://mrsrl.stanford.edu/~brian/bloch/

Ry = [cos(phi) 0 sin(phi);0 1 0;-sin(phi) 0 cos(phi)];

