function coVarMat = CalculatePFCovariance(peakData, threshold, omeStep, instr, options)
% CalculatePFCovariance:  This function calculates the covariance matrix of
% the pole figure.  Each gS is weighted by the corresponding intensity.
% 
% USAGE:
%   coVarMat = CalculatePFCovariance(peakData, threshold, omeStep, instr, options)
% 
% AUTHOR: Timothy Long
% 
% INPUTS:
%   peakData is a 1 x 1 structure:
%       .tth0 is 1 x 1:
%           The unstrained two theta value of the peak
%       .eta0 is 1 x 1:
%           The eta center of the peak as determined by peak fitting
%       .ome0 is 1 x 1:
%           The omega center of the peak as determined by peak fitting
%       .IList is n (number of integration points) x 1:
%           The integrated (over tth) intensity of each of the 'n' points
%       .etaList is n x 1:
%           The eta value of each of the 'n' points
%       .omeList is n x 1:
%           The omega value of each of the 'n' points
% 
%   threshold is 1 x 1:
%       A threshold used to convert the pole figure to black and white to
%       identify the one 'spot' assosiated with one diffraction peak from
%       one grain.
%
%   omeStep is 1 x 1:
%       The size of each omega step, i.e. the difference in omega between
%       two adjacent detector frames.
%
%
%   options is a 1 x 1 structure:
%       The options structure.  Currently not used.
%
%   NOTES:
%       Started 2015/1/14
%
%       Based on the work of Robert Carson.


% unpack the structure
fitTth = peakData.tth0;
fitEta = peakData.eta0;
fitOme = peakData.ome0;


% convert pfData to an image in omega-eta space
[ome,index1] = sort(peakData.omeList);
tth = peakData.tthList(index1);
eta = peakData.etaList(index1);
I = peakData.IList(index1);
minOme = min(ome);
maxOme = max(ome);
omeRange = maxOme-minOme;
width = uint16(omeRange / omeStep + 1);
height = uint16(length(ome)/width);

omeMesh = reshape(ome,height,width);
tthMesh = reshape(tth,height,width);
etaMesh = reshape(eta,height,width);
IMesh = reshape(I,height,width);

% find the spot closest to the peak fit eta and ome
% spotIndex = IDSpot(IMesh, etaMesh, omeMesh, threshold, nominalCoords)
index2 = IDSpot(IMesh, etaMesh, omeMesh, threshold, [fitEta,fitOme]);

spotI = IMesh(index2);
spotTth = tthMesh(index2);
spotEta = etaMesh(index2);
spotOme = omeMesh(index2);

gS = qFromTthEtaOme(instr,[fitTth;spotTth],[fitEta;spotEta],[fitOme;spotOme]);
gS0 = gS(1,:);
gSspot = gS(2:end,:);

% center the data so the average is 0
deltaGS = [gSspot(:,1)-gS0(1), gSspot(:,2)-gS0(2), gSspot(:,3)-gS0(3)];

% calculate the covariance matrix of gS
coVarMat = zeros(3);
totalI = sum(spotI);
for ii = 1:length(spotI)
    x = (deltaGS(ii,:)'*deltaGS(ii,:));
    coVarMat = coVarMat + (spotI(ii)/totalI * x);
end

end