function pfData = MainSingleGrainPFProcessing(images, instr, material, cakeParms, grainData, omeBnds, omeStep, options )
% MainSingleGrainPFProcessing: This function builds single crystal pole
% figures for peaks from a single grain in an aggregate of grains.  This
% function requires the data to have already been run through
% MainGrainFinder and MainGrainProcessing.
% 
% Author: Timothy Long (probably also Darran Pagan and Mark Obstalecki)
% 
% USAGE: pfData = MainSingleGrainPFProcessing(images, instr, material, ...
%                 cakeParms, grainData, omeBnds, omeStep, options )
% 
% INPUTS:
%   images is instr.det.pixelDim(1) x instr.det.pixelDim(2) x n
%       A stack of all the images from the detector.  WARNING: This stack
%       is very large in RAM.
% 
%   instr is a 1x1 structure
%       A structure containing parameters describing the physical
%       intsurmentation used in the experiment.  Output by APEX_RD_Start
%       and modified by MainCalibrantProcessing.
% 
%   material is a 1x1 structure
%       A structure containing the parameters describing the crystal
%       structure of the sample material.  Output by APEX_RD_Start
%
%   cakeParms is a 1x1 structure
%       A structure containing the parameters defining the polar
%       integration (cake) of the detector images.  Output by APEX_RD_Start
%       and modified by MainCalibrantProcessing.
% 
%   grainData is a 1 x m structure
%       A structure containing the orientation, position, strain tensors,
%       and fitData for m grains.  Output by MainGrainProcessing
% 
%   omeBnds is 1x2
%       The lower and upper bounds for omega in the image stack
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
%       .hklList is k x 3:
%           A list of the hkl's to build pole figures for, if those hkl's
%           diffract for the input grains.  Default is material.hkls.
%
%       .dome is 1x1
%           The number of frames on each side of the center frame to use
%
% OUTPUTS:
%   pf is a 1 x m structure  
%       Contains the information for pole figures for each of the m grains.
%       Each grain can have pole figures for multiple peaks depending on
%       the input options.
%           .peaks
%
% NOTES:
%   Started 2015_June_28
%
%   To Darren, Mark, and/or Chris:  If you're editing this and want to add
%   your name to the author list (and I haven't yet done so) feel free to
%   add yourself.


% Set default options

hklList = material.hkls;
dome = 4;

% Overwrite default options if pass in
if isfield(options,'hklList')
    hklList = options.hklList;
end

if isfield(options,'dome')
    dome = options.dome;
end

numGrains = size(grainData,2);

numPeaks = size(hklList,1);

numImgs=size(images,3);

% Enable multiprocessing 
% if(isempty(gcp('nocreate'))) %checking to see if my pool is already open
%     parpool
% end

pfData = struct();

% parfor ii = 1:numGrains
for ii = 1:numGrains
    for jj = 1:numPeaks
        
        % check if peak is predicted to diffract onto the detector in
        % omeBnds
        [~,index1] = ismember(hklList(jj,:),grainData(ii).fitData.hkl,'rows');
        
        % create temporary storage variable
        peak = struct();
        peak.hkl = hklList(jj,:);
        
        %%% DEBUG %%%
        disp(peak.hkl)
        disp(index1)
        disp(grainData(ii).fitData.frameno(index1))
        disp(grainData(ii).fitData.ome(index1))
        %%%%%%%%%%%%%
        
        % check that the peak is a good one
        check1 = ~isempty(index1);
        check2 = ismember(index1,grainData(ii).fitData.pksToUse);
        
        % I'm not passing in the whole image stack, so avoid peaks near the
        % edges
        check3 = (omeBnds(1)+omeStep*(dome+1)) < grainData(ii).fitData.ome(index1);
        
        check4 = grainData(ii).fitData.ome(index1) < (omeBnds(2)-omeStep*dome);
        
        %%% DEBUG %%%
        disp(check1 && check2 && check3 && check4)
        disp('--------')
        disp(' ')
        %%%%%%%%%%%%%
        
        if(check1 && check2 && check3 && check4)
            % get the center of the peak
            posDet = [grainData(ii).fitData.pixCenter(index1,:),grainData(ii).fitData.frameno(index1)];
            
            
            % entryNo = 1; % This variable is only used if preallocation is used.
            
            
            % for each frame, radially integrate the peak
            for kk = posDet(3)-dome:posDet(3)+dome
                
                %%% debug %%%
                
                disp(kk)
                
                %%%%%%%%%%%%%
                
                % integrate peak
                integratedPeak = IntegratePeak(images(:,:,kk),posDet(1:2),instr,cakeParms,options);

                % find the number of data points
                numPoints = length(integratedPeak.etaList);
                
                % find rought omega for the frame
                omeR=FindOmeRough(omeBnds(1),omeBnds(2),omeStep,kk);
                
                % assign the rough omega value to all eta points
                omeList = omeR*ones(numPoints,1);
                
                % assign the peak's center two theta value to all eta points
                tthList = grainData(ii).fitData.tth(index1)*ones(numPoints,1);
                
                %%% Note:
                %   The following is probably rather slow and needs to be
                %   changed to use preallocation (like MainSCPFProcessing)
                %
                
                if(kk == posDet(3)-dome)
                    % create fields on the first pass
                    
                    % calculate and save the sample scattering vectors
                    peak.gS = qFromTthEtaOme(instr,tthList,...
                    integratedPeak.etaList,omeList);
                
                    % save the radially integrated intensity for those
                    % scattering vectors
                    peak.I = integratedPeak.etaIList;
                    
                else
                    % concatenate to fields on following passes
                    
                    % calculate and save the sample scattering vectors
                    peak.gS = cat(1,peak.gS,qFromTthEtaOme(instr,tthList,...
                        integratedPeak.etaList,omeList));

                    % save the radially integrated intensity for those
                    % scattering vectors
                    peak.I = cat(1,peak.I,integratedPeak.etaIList);
                end
                
                %%%
                
                % This section is only needed if preallocation of gS and I
                % is used
                
                % increment the index by the number of data points
                % entryNo = entryNo + numPoints;
            end
            
            % write data to output variable
            pfData(ii).peak(jj) = peak;
        end
    end
end


end