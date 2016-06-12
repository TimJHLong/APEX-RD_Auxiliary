function q = GenerateQuatGridLR(n,ii)
% GenerateQuatGridLR: A function to generates points in an array over all
% of quaternion orientation space one point at a time.  This should make it
% compatable with parfor loops.
% 
% USAGE: q = GenerateQuatGridLR(n,ii)
% 
% AUTHOR: Timothy Long
% 
% INPUTS:
%   n is a scalar:
%       The number of points per edge of the unit tesseract.  A larger
%       value for n produces a denser array of points.
% 
%   ii is a scalar:
%       The number of the point in the grid to return.
%
% OUTPUT:
%   q is 1 x 4:
%       A unit quaternion.
%
% NOTES:
%   Started 2016/May/30
%
%   Needs testing before being uses.


dist = 2/n;
list = (-1:dist:1);

rowNum = modmin1(ii,(n+1));
colNum = modmin1(ceil(ii/(n+1)), (n+1));
sheetNum = modmin1(ceil(ii/(n+1)^2),(n+1));
faceNum = ceil(ii/(n+1)^3);

% disp(['Row #' num2str(rowNum)])

switch faceNum
    case 1 % x = 1 face
        q = [1, list(rowNum), list(colNum), list(sheetNum)];
    case 2 % x = -1 face
        q = [-1, list(rowNum), list(colNum), list(sheetNum)];

    case 3 % y = 1 face
        q = [list(rowNum), 1, list(colNum), list(sheetNum)];
    case 4 % y = -1 face
        q = [list(rowNum), -1, list(colNum), list(sheetNum)];

    case 5 % z = 1 face
        q = [list(rowNum), list(colNum), 1, list(sheetNum)];
    case 6 % z = -1 face
        q = [list(rowNum), list(colNum), -1, list(sheetNum)];

    case 7 % w = 1 face
        q = [list(rowNum), list(colNum), list(sheetNum), 1];
    case 8 % w = -1 face
        q = [list(rowNum), list(colNum), list(sheetNum), -1];
end

% project the point onto the unit 4-sphere
q = q/norm(q);

end


function m = modmin1(n,base)
    % a helper function that is like mod, but has an output range of
    % [1,base] instead of [0,base-1].  Useful for dealing with one based
    % indexing.
    m = mod(n,base);

    if m ==0
        m = base;
    end
end