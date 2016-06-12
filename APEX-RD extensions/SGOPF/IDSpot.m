function spotIndex = IDSpot(IMesh, etaMesh, omeMesh, threshold, nominalCoords)
% IDSpot:  This function takes a single color image containing one or
% more 'hot spots' with intensities greater than the threshold value and
% finds the one spot with a center of intensity closest to the nominal x
% and y coordinates.
% 
% USAGE: spotIndex = IDSpot( image, threshold, nominalCoords )
% 
% AUTHOR: Timothy Long
%
% INPUTS:
%   image is n (y-dimension) x m (x-dimension):
%       A 2D matrix representing a single color image or intensities
% 
%   threshold is 1 x 1:
%       A value used to convert the image to a logical black and white
%       image.  Values below the threshold become 0, while values equal to
%       or greater than the threshold become 1.
% 
%   nominalCoords is 2 x 2:
%       The approximate x & y coordinates, in pixels within the image, of
%       the center of the spot to pick out.
% 
% OUTPUTS:
%   spotIndex is p x 1:
%       The linear indicies of the pixels in image that fall within the
%       spot.
% 
% NOTES:
%   Started 2016/1/14
%   
%   The first axis (vertical when the matrix is displayed) is 'y' and the
%   second axis (horizontal when the matrix is displayed) is 'x'.
%
%   This code is intended to be used to pick out the spot assosiated with a
%   particular peak from one hkl from one grain in a multigrain pole figure


bwImage = IMesh >= threshold;

% find all the spots in the image
spots = bwlabel(bwImage);

% calculate center of intensities for the spots
spotsData = struct();
for ii = 1:max(max(spots))
    index1 = (spots==ii);
    tempI = IMesh(index1);
    tempEta = etaMesh(index1);
    tempOme = omeMesh(index1);
    
    totalI = sum(tempI);
    etaCenter = sum(tempI.*tempEta)/(totalI);
    omeCenter = sum(tempI.*tempOme)/(totalI);
    
    dist = ((nominalCoords(1)-etaCenter)^2 + (nominalCoords(2)-omeCenter)^2)^0.5;
    
    spotsData(ii).etaCenter = etaCenter;
    spotsData(ii).omeCenter = omeCenter;
    spotsData(ii).dist = dist;
    
    
    %%% DEBUG %%%
%     disp(['Spot #' num2str(ii)])
%     disp(['Distance = ' num2str(spotsData(ii).dist)])
%     tempImg = zeros(size(IMesh));
%     tempImg(spots==ii) = IMesh(spots==ii);
%     imagesc(tempImg)
%     keyboard
    %%%%%%%%%%%%%
end

% find the spot closest to the center
distList = [spotsData.dist];
[~,index2] = min(distList);

spotIndex = find(spots == index2);

end