function grainList = findRepeatedGrains(grainData,distTol,angleTol)
% THIS FUNCTION IS ODFPF DEPENDENT!!
%
% findRepeatedGrains: This function finds grains that are within some
% physical distance and an angular distance in orientation space of each
% other and removes all but the one that has the largest number of well fit
% peaks. This is a first pass at finding repeated hits on the same grain
% and is not very sophisticated.
%
% AUTHOR: Timothy Long
%
% USAGE: grainList = findRepeatedGrains(grainData,distTol,angleTol)
%
% INPUTS:
%   grainData is a 1 x n structure
%       The output from fine peak processing (MainPeakProcessing or its low
%       RAM version)
% 
%   distTol is 1x1 scalar
%       
% 
% 





warning('This function is ODFPF dependent!')


end