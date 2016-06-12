function [trajStart, trajEnd, trajMag] = AutoCalcPFTrajectory(peakData, threshold, fitEta, fitOme, omeStep, options)
% autoCalcPFTrajectory:  This function converts the pole figure to an
% image, uses the image processing toolbox to find all continuous bright
% spots in the image, finds the one with its center of intensity closest to
% the fit peak center, and then finds the longest straight line that can
% fit in the peak.
% 
% USAGE:
%   [trajStart, trajDir, trajMag] = AutoCalcPFTrajectory(peakData,...
%       threshold, fitEta, fitOme, omeStep, options)
% 
% AUTHOR: Timothy Long (maybe Darren Pagan)
% 
% INPUTS:
%   peakData: a  1 x 1 structure with the fields:
%       .IList is n x 1:
%           The intensity of each point on the pole figure.
%       .etaList is n x 1:
%           The eta values (in degrees) of each point
%       .omeList is n x 1:
%           The omega values (in degrees) of each point
% 
%   threshold is 1 x 1:
%       A threshold used to convert the pole figure to a binary image for
%       processing.  Points with intensities below the threshold are set to
%       0, and points equal to or above it are set to 1.  The value is a
%       percentage relative to the intensity of the maxima closest to the
%       fit peak center.
% 
%   fitEta is 1 x 1:
%       The eta position of the peak measured by peak fitting.
% 
%   fitOme is 1 x 1
%       The eta position of the peak measured by peak fitting.
% 
%   omeStep is 1 x 1:
%       The size of each omega step, i.e. the difference in omega between
%       two adjacent detector frames.
%
%   etaStep is 1 x 1:
%       The difference in eta between any two integration points with the
%       same omega value.
%
%   options is a 1 x 1 structure:
%       The options structure.  Currently not used.
% 
% OUTPUTS:
%   trajStart is 2 x 1:
%       The eta and omega positions of the first point on the trajectory
%
%   trajEnd is 2 x 1:
%       The eta and omega positions of the last point on the trajectory
%
%   trajMag is 1 x 1:
%       The magnitude of the trajectory vector in degrees (eta & omega
%       coordinates).
% 
% NOTES:
%   Started 2016/1/14
%
%   Unlike MeasurePFTrajectory, this function uses omega and eta
%   coordinates instead of gS_x,gS_y,gS_z.
%
% ACKNOWLEDGEMENTS:
%   This function is based on the work of Darren Pagan.

% defaults
debug = false;
numTrajs = 1;
trajAngle = 15;
useIntI = false;

if exist('options','var')
    if isfield(options,'numTrajs')
        numTrajs = options.numTrajs;
    end
    
    if isfield(options,'trajAngle')
        trajAngle = options.trajAngle;
    end
    
    if isfield(options,'debug')
        debug = options.debug;
    end
    
    if isfield(options,'useIntI')
        useIntI = options.useIntI;
    end
end

% convert pfData to an image in omega-eta space
[ome,index1] = sort(peakData.omeList);
eta = peakData.etaList(index1);
I = peakData.IList(index1);

minOme = min(ome);
maxOme = max(ome);
omeRange = maxOme-minOme;
width = uint16(omeRange / omeStep + 1);
height = uint16(length(ome)/width);

omeMesh = reshape(ome,height,width);
etaMesh = reshape(eta,height,width);
IMesh = reshape(I,height,width);

% find the maximum closest to the fit point and scale the entire image by
% that intensity before thresholding
index4 = findLocalMax(IMesh, 5, 5);
etaMaxes = etaMesh(index4);
omeMaxes = omeMesh(index4);
IMaxes = IMesh(index4);
deltas = [etaMaxes-fitEta, omeMaxes-fitOme];
dist = sum(deltas.^2,2).^0.5;
[~,index5] = min(dist);
IMesh = IMesh./IMaxes(index5) * 100;

% find the spot closest to the peak fit eta and ome
center = [fitEta, fitOme];
index2 = IDSpot(IMesh, etaMesh, omeMesh, threshold, center);

spotI = IMesh(index2);
spotEta = etaMesh(index2);
spotOme = omeMesh(index2);

%%%
%
% EXPERIMENTAL: find highest total intensity line in the PF
%
if useIntI
trajStart = zeros(2,numTrajs);
trajEnd = zeros(2,numTrajs);
trajMag = zeros(numTrajs,1);
trajIntI = zeros(numTrajs,1);

[~,index3] = max(spotI);
startPt = [spotEta(index3);spotOme(index3)];

for nt = 1:numTrajs
    % avoid endPt = startPt
    for ii = [(1:index3-1), (index3+1:length(spotI))]
        endPt = [spotEta(ii);spotOme(ii)];
        vec = endPt - startPt;
        mag = norm(vec);
        numEl = ceil(mag/omeStep);
        lE = mag/numEl;
        pts = [(startPt(1) + vec(1)/numEl*(0:numEl));
               (startPt(2) + vec(2)/numEl*(0:numEl))];
        
        % interpolate the intensity for each point on the line
        IE = zeros(numEl+1,1);
        IE(1) = spotI(index3);
        IE(numEl+1) = spotI(ii);
        for jj = 2:numEl
            % find the 'element' the point lies in
            lInd1 = etaMesh <= pts(1,jj);
            lInd2 = omeMesh <= pts(2,jj);
            lInd3 = lInd1 & lInd2;
            
            [ind1,ind2] = find(lInd3);
            ind3 = max(ind1);
            ind4 = max(ind2);
                        
            checkSz = size(IMesh);
            if ind3>checkSz(1) || ind3<=0 || ind4>checkSz(2) || ind4<=0
                disp('Error, point in line falls on edge of pole figure')
                disp(['Line #' num2str(ii)])
                disp(['El #' num2str(jj)])
                disp(ind3)
                disp(ind4)
                keyboard
            end
            
            % if a line runs along the edge, use the edge elements instead
            % of looking for an element that doesn't exist.
            if ind3==checkSz(1);
                ind3 = ind3-1;
            elseif ind4==checkSz(2);
                ind4 = ind4-1;
            end
            
            % use bilinear interpolation to estimate the intensity at the
            % point.
            etaStep = etaMesh(ind3+1,ind4)-etaMesh(ind3,ind4);
            x = [(pts(1)-etaMesh(ind3,ind4))/etaStep;
                 (pts(2)-omeMesh(ind3,ind4))/omeStep];
            
            IE(jj) = ((1-x(2))*(1-x(1)) * IMesh(ind3,ind4) + ...
                     (1-x(2))*x(1) * IMesh(ind3+1,ind4) + ...
                     x(2)*(1-x(1)) * IMesh(ind3,ind4+1) + ...
                     x(2)*x(1) * IMesh(ind3+1,ind4+1));
        end
        
        % apply trapezoid rule to integrate intensity
        weight = 2*ones(size(IE));
        weight(1) = 1;
        weight(end) = 1;
        intI = weight' * IE * lE;
        
        % check angle with other trajectories
        if nt > 1
            tempVecs = (trajEnd(:,1:nt-1)-trajStart(:,1:nt-1));
            dots = vec' * tempVecs./(norm(vec)*sum(tempVecs.^2,1).^0.5);
            angles = 180/pi * acos(dots);
        else
            angles = 180;
        end
        
        % record the trajectory
        if (intI > trajIntI(nt)) && (sum(angles>trajAngle)==length(angles))
            trajStart(:,nt) = startPt;
            trajEnd(:,nt) = endPt;
            trajMag(nt) = mag;
            trajIntI(nt) = intI;
        end
    end
end
else

%%% 
% 
% find the longest vector from the fit center to a point in the spot
%
trajStart = zeros(2,numTrajs);
trajEnd = zeros(2,numTrajs);
trajMag = zeros(numTrajs,1);

% test all possible vectors and store the largest
[~,index3] = max(spotI);
for nt = 1:numTrajs
    for ii = 1:length(spotI)
        startPt = [spotEta(index3);spotOme(index3)];
        endPt = [spotEta(ii);spotOme(ii)];

        vec = endPt - startPt;
        mag = norm(vec);
        
        if nt > 1
            tempVecs = (trajStart(:,1:nt-1)-trajEnd(:,1:nt-1));
            dots = vec' * tempVecs./(norm(vec)*sum(tempVecs.^2,1).^0.5);
            angles = 180/pi * acos(dots);
        else
            angles = 180;
        end
        
        if(mag > trajMag(nt)) && (sum(angles>trajAngle)==length(angles))
            trajStart(:,nt) = startPt;
            trajEnd(:,nt) = endPt;
            trajMag(nt) = mag;
        end
    end
end
end

end