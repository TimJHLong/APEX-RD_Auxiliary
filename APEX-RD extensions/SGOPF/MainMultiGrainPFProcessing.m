function pfData = MainMultiGrainPFProcessing(images, instr, cakeParms, grainData, omeBnds, omeStep, options)
% MainMultiGrainPFProcessing: This function processes the input images to
% produce the data for pole figures for all grains in grainData.  This code
% does not currently support split omega bounds, ie front & back sweeps
% with a missed area due to the loadframe's posts.
%
% USAGE:
%   pfData = MainMultiGrainPFProcessing(images, instr, cakeParms,
%   grainData, omeBnds, omeStep, options)
%
% AUTHOR: Timothy Long & Darren Pagan (used code from MainSCPFProcessing)
%
% INPUTS:
%   images is instr.det.pixelDim(1) x instr.det.pixelDim(2) x r (# of images)
%       A stack of all the detector images in one load sweep
%
%   instr is a 1x1 structure
%       A structure containing parameters describing the physical
%       intsurmentation used in the experiment.  Output by APEX_RD_Start
%       and modified by MainCalibrantProcessing.
%
%   cakeParms is a 1x1 structure
%       A structure containing the parameters defining the polar
%       integration (cake) of the detector images.  Output by APEX_RD_Start
%       and modified by MainCalibrantProcessing.
%           .peakParms is a 1x1 structure
%               WARNING! Changing cakeParms.peakPars will impact the
%               quality of peak fits done by MainGrainProcessing!  Remember
%               to reset the values before attempting to fit peaks.
%                   .radIntLength is 1x1
%                       The radial distance, in pixels, over which the
%                       images are integrated.  Increases this to include
%                       a larger two theta range in each point.
%                   .etaIntLength is 1x1
%                       The azimuthal distance, in degrees, over which the
%                       images are radially integrated.  Increase this for
%                       a "taller" pole figure.
%
%   grainData is a 1 x p structure
%       A structure containing the orientation, position, strain tensors,
%       and fitData for p grains.  Output by MainGrainProcessing
% 
%   omeBnds is m x 2
%       The lower and upper bounds for omega for each of the m sweeps.
%       Note that this is the sames as for MainGrainProcessingLowRam.  Most
%       APEX_RD functions will accept this format, but this function will
%       not accept the global limits like most other functions will.
% 
%   omeStep is 1x1
%       The rotation angle in omega over which each image is captured.
%       Also the rotation angle in omega between images.
% 
%   options is a 1x1 structure
%       Contains options that control how the code executes.
%       .matchRadIntScaling is 1x1:
%           If this options is enable, the value in
%           cakeParms.peakParms.etaIntLength will be ignored and the value
%           from cakeParms.peakParms.radIntLength will be used instead.
%           Default is 0.
%
%       .dome is 1x1
%           The number of frames on each side of the center frame to use
%
% OUTPUTS:
%   pfData is a 1 x numGrains structure
%       .peaks is a 1 x numPeaks structure
%           .etaIList is n (# of integrated points) x 1
%               A list of the radially integrated intensities for a peak
%           .etaList is n x 1
%               A list of the eta values at which the peak was integrated
%           .gS is n x 3
%               The scattering vector for each integrated point
%           .hkl is 1 x 3
%               The Miller indicies of the peak 
% 
% NOTES:
%   Started 2015_June_30
%
%   Requires all the fastsweeps be loaded into RAM at once
%
%   Does not currently support split omega bounds, bounds must be
%   continuous.
%
%   Darren Pagan is listed as an author to credit Tim's use of his code in
%   this function.  However, Darren has not seen this function and is not
%   responsible for it quality or lack thereof.



% default options
dome = 4;
filterPeaks = false;

% overwrite options if passed in
if(isfield(options,'dome'))
    dome = options.dome(1);
end

if(isfield(options,'filterPeaks'))
    filterPeaks = options.filterPeaks;
end


% find parameters
numGrains = size(grainData,2);


% initialize storage variables
pfData = struct();


for kk = 1:numGrains
            
    disp(['Processing grain #' num2str(kk)])

    omeLow = omeBnds(1) + (dome+1)*omeStep;
    omeHigh = omeBnds(2) - (dome+1)*omeStep; 

    % find all the peaks inside the omega bounds
    index1 = grainData(kk).fitData.ome > omeLow;
    index2 = grainData(kk).fitData.ome < omeHigh;
    index3 = index1 & index2;
    
    if filterPeaks
        % convert pksToUse to a logical array
        index4 = zeros(size(grainData(kk).fitData.ome,1),1);
        index4(grainData(kk).fitData.pksToUse) = 1;
        index4 = logical(index4);


        % select only the well fit peaks that are fully within the bounds
        frameno = grainData(kk).fitData.frameno(index3 & index4);
        pixCoords = grainData(kk).fitData.pixCenter(index3 & index4,:);
        tths = grainData(kk).fitData.tth(index3 & index4);
        hkls = grainData(kk).fitData.hkl(index3 & index4,:);
    else
        
        % select only the peaks that are fully within the bounds
        frameno = FindFrameNo(grainData(kk).fitData.ome(index3),omeBnds,omeStep);
        pixCoords = grainData(kk).fitData.pixCenter(index3,:);
        tths = grainData(kk).fitData.tth(index3);
        hkls = grainData(kk).fitData.hkl(index3,:);
    end
    
    % integrate each peak
    for ll = 1:size(hkls,1)
        disp(['Processing peak number ' num2str(ll)])
        % create local storage variable and preallocate sizes
        peak = struct();
        numPoints = 100 * (2*dome+1) * cakeParms.peakParms.etaIntLength * cakeParms.peakParms.binscaling;
        peak.IList = zeros(numPoints,1);
        peak.etaList = zeros(numPoints,1);
        peak.omeList = zeros(numPoints,1);
        peak.tthList = zeros(numPoints,1);
        peak.gS = zeros(numPoints,3);
        peak.hkl = hkls(ll,:);


        %%% DEBUG %%%
        if filterPeaks
            check1 = index3 & index4;
        else
            check1 = index3;
        end
        debugLogical1 = ismember(peak.hkl,grainData(kk).fitData.hkl(check1,:));
        if(~debugLogical1)
            disp([kk, ll])
            disp(peak.hkl)
            error('Selected HKL does not exist.')
        end
        %%%%%%%%%%%%%

        % create a small stack of images so the integration can be
        % parallelized
        tempImages = images(:,:,frameno(ll)-dome:frameno(ll)+dome);
        tempFramenos = frameno(ll)-dome:frameno(ll)+dome;


        % integrate each frame
        for mm = 1:(2*dome+1)
            integratedPeak = IntegratePeak(tempImages(:,:,mm),pixCoords(ll,:),instr,cakeParms,options);


            % build lists to calculate gS with
            etaList = integratedPeak.etaList;


            % assign the unstrained tth to all points in the peak
            tthList = ones(length(etaList),1)*tths(ll);


            % finding omeR
            omeR = FindOmeRough(omeBnds(1),omeBnds(2),omeStep,tempFramenos(mm));


            % find scattering vector
            omeList = ones(length(etaList),1)*omeR;
            qList = qFromTthEtaOme(instr,tthList,integratedPeak.etaList,omeList);


            % write to output variable
            peak.IList = [peak.IList; integratedPeak.etaIList];
            peak.tthList = [peak.tthList; tthList];
            peak.etaList = [peak.etaList; etaList];
            peak.omeList = [peak.omeList; omeList];
            peak.gS = [peak.gS; qList];
        end


        % save the data to the output structure
        pfData(kk).peaks(ll) = peak;

    end
end


end