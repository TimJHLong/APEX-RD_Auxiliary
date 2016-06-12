function [outliers, average, StD] = findOutliers(data,threshold, numSTD, options)
% findOutliers:  This function takes a set of scalar data in the form of a
% vector and finds outliers based on an absolute threshold and a number of
% standard deviations away from average.  
%
% USAGE: outliers = findOutliers(data,threshold, numSTD, options)
%
% AUTHOR: Timothy Long 
%
% INPUTS:
%   data is n x 1
%       A vector containing n scalar values to find the outliers in
% 
%   threshold is 1 x 1
%       Entries in 'data' that have absolute values larger than 'threshold'
%       will be marked as outliers regardless of the average and standard
%       deviation of the data values.  This values should always be
%       positive.  If no threshold is desired, input NaN or an empty
%       matrix.
% 
%   numSTD is 1 x 1
%       The number of standard deviations away from average a values is
%       allowed to be without being flagged as an outlier.  Values such
%       that |value - average| > numSTD * (standard deviation) are flagged
%       as outliers.
% 
%   options is a 1 x 1 structure
%       This structure contains assorted options that control how the code
%       executes.  Options used by this function are:
%           .excludeCurrentPoint is 1 x 1
%               If set to 1, the point being tested will be excluded from
%               the calculation of the average and standard deviation.
%               This slows down execution, but tightens tolerances.  If set
%               to 0, the average and std of all points will be used. 
%               Default is 0.
%           .numPasses is 1 x 1
%               The number of times to repeat the test for outliers on the
%               data set, each time excluding all points marked as outliers
%               by all previous loops.  This value must be positive.
%               Increasing the number of passes helps filter out smaller
%               outliers in a data set that also contains one or two large
%               outliers.  Default is 1. 
%
% OUTPUTS:
%   outliers is n x 1
%       A vector of "logical" values (0s and 1s).  The i'th entry is 1 if
%       the i'th entry of data is found to be an outlier.
%
%   avgerage is 1 x 1
%       The average value of all non-outlier points
%
%   StD is 1 x 1
%       The standard deviation of all non-outlier points
%
% Notes:
%   Started 2015_May_7


% set default options
numPasses = 1;
exCP = 0;

% overwrite options if passed in
if isfield(options,'excludeCurrentPoint')
    exCP = options.excludeCurrentPoint;
end

if isfield(options,'numPasses')
    numPasses = options.numPasses;
end


% initialize storage variable
n = length(data);
outliers = zeros(n,1);


% if threshold is a number, find any values in abs(data) greater than
% threshold and flag them as outliers.
if ~isempty(threshold) && ~isnan(threshold)
    % note the 'greater than,' not 'greater than or equal to'
    outliers = abs(data)>threshold;
end


% filter by standard deviation 
for ii=1:numPasses
    
    if exCP % exclude the test point from the average and STD calculations
        
        % itterate over all points
        for jj = 1:n
            % re-test previous outliers, just to be sure
            
            % create a copy of the data that excludes previous outliers and
            % the point being tested
            tempOut = outliers;
            tempOut(jj) = 1;
            tempData = data(~tempOut);

            % calculate the average and std of the included data
            len = length(tempData);
            avg = sum(tempData)/len;
            std = (sum((tempData - avg).^2)/(len-1))^0.5;

            % check if the test point is an outlier
            if (data(jj)-avg) > numSTD * std
                outliers(jj) = 1;
            end
                

            
        end
        
    else % use all data point to calculate average and STD
        
        % filter out already flagged outliers
        tempData = data(~outliers);
        len = length(tempData);
        
        % calculate average and standard deviation
        avg = sum(tempData)/len;
        std = ( sum((tempData-avg).^2) / (len-1) )^0.5;
        
        % look for values more than numSTD*std away from the average
        outsideSTD = (data-avg) > numSTD * std;
        
        % combine lits ot outliers with existing list.  The use of array or
        % ( | )  is to keep any previous outliers flagged as such.
        outliers = outliers | outsideSTD;
    end
    
end

% calculate final average and standard deviation
fData = data(~outliers); % filtered data
average = sum(fData)/length(fData);
StD = ( sum((fData - average).^2)/(length(fData)-1) )^0.5;

end