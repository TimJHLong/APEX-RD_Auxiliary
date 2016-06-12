function s = cubicSymmetry
% cubicSymmetry:  This function return the 24 cubic symmetry rotation
% matricies.
%
% USAGE: s = cubicSymmetry;
%
% INPUTS: NONE
%
% OUTPUTS:
%   s: s is a 3 x 3 x 24 matrix where each 3 x 3 x 1 slice is one of the 24
%   cubic symetries represented as a rotation matrix.
%
% AUTHOR: Timothy Long
%
% Notes:
%   Started 2015 Feb 23
%
%   Source: Engler, O. & V. Randle.  "Texture Analysis: Macrotexture,
%   Microtexture & Orientation Mapping".  Gordon and Breach Science
%   Publishers (2000).


s = zeros(3,3,24);

%24 cubic sysmetries
s(:,:,1) = [1 0 0; 0 1 0; 0 0 1];
s(:,:,2) = [0 0 -1; 0 -1 0; -1 0 0];
s(:,:,3) = [0 0 -1; 0 1 0; 1 0 0];
s(:,:,4) = [-1 0 0; 0 1 0; 0 0 -1];
s(:,:,5) = [0 0 1; 0 1 0; -1 0 0];
s(:,:,6) = [1 0 0; 0 0 -1; 0 1 0];
s(:,:,7) = [1 0 0; 0 -1 0; 0 0 -1];
s(:,:,8) = [1 0 0; 0 0 1; 0 -1 0];
s(:,:,9) = [0 -1 0; 1 0 0; 0 0 1];
s(:,:,10) = [-1 0 0; 0 -1 0; 0 0 1];
s(:,:,11) = [0 1 0; -1 0 0; 0 0 1];
s(:,:,12) = [0 0 1; 1 0 0; 0 1 0];
s(:,:,13) = [0 1 0; 0 0 1; 1 0 0];
s(:,:,14) = [0 0 -1; -1 0 0; 0 1 0];
s(:,:,15) = [0 -1 0; 0 0 1; -1 0 0];
s(:,:,16) = [0 1 0; 0 0 -1; -1 0 0];
s(:,:,17) = [0 0 -1; 1 0 0; 0 -1 0];
s(:,:,18) = [0 0 1; -1 0 0; 0 -1 0];
s(:,:,19) = [0 -1 0; 0 0 -1; 1 0 0];
s(:,:,20) = [0 1 0; 1 0 0; 0 0 -1];
s(:,:,21) = [-1 0 0; 0 0 1; 0 1 0];
s(:,:,22) = [0 0 1; 0 -1 0; 1 0 0];
s(:,:,23) = [0 -1 0; -1 0 0; 0 0 -1];
s(:,:,24) = [-1 0 0; 0 0 -1; 0 -1 0];

end