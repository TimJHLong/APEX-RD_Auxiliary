function pfData = MainMultiGrainPFProcessingLR(imgInfo, instr, cakeParms, allLoadData, omeBnds, omeStep, options)
% MainMultiGrainPFProcessingLR: This function processes the images defined
% by imgInfo to produce the data for pole figures for all grains in 
% grainData across all load steps.  This code loads one fastsweep at a time
% to avoid using virtual memory, which is very slow.
% 
% USAGE: pfData = MainMultiGrainPFProcessingLR(imgInfo, instr, ...
%           cakeParms, allLoadData, omeBnds, omeStep, options)
% 
% AUTHOR: Timothy Long (probably also Darran Pagan and Mark Obstalecki)
% 
% INPUTS:
%   imgInfo is a structure containing the information needed to load the
%   images.  Currently only supports .ge2 image files (ex from CHESS F2)
%       .dir is a string 
%           complete path to the fast sweep file, e.g. 'C:/data/'  All load
%           steps' images must be in the same directory (for now)
%
%       .stem is a string
%           stem of the files.  This must be the same across all files.
%           For CHESS farfield data, it should be 'ff'.
%
%       .sweepNumArray is n x m
%           matrix where each row gives the m image numbers for one of the
%           n loadsteps.  Must be square, like all matlab matricies.
%
%       .numDark is 1 x 1
%           number of dark images at the front of the fast sweep files     
% 
%       .numPerStack is 1 x 1
%           numbers of data images in the .ge2 file
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
%   allLoadData is a 1 x n (# of load steps) structure
%       Contains the grain data for all load steps.  Has fields:
%       
%       .grainData is a 1 x p structure
%           A structure containing the orientation, position, strain
%           tensors, and fitData for p grains.
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
%       .hklList is k x 3:
%           A list of the hkl's to build pole figures for, if those hkl's
%           diffract for the input grains.  Default is material.hkls.
%
%       .dome is 1x1
%           The number of frames on each side of the center frame to use
%
% OUTPUTS:
%   pfData is a 1 x numLoadSteps structure
%       .grains is a 1 x numGrains structure
%           .peaks is a 1 x numPeaks structure
%               .etaIList is n (# of integrated points) x 1
%                   A list of the radially integrated intensities for a peak
%               .etaList is n x 1
%                   A list of the eta values at which the peak was integrated
%               .gS is n x 3
%                   The scattering vector for each integrated point
%               .hkl is 1 x 3
%                   The Miller indicies of the peak 




%
if(~strcmp(instr.det.ext,'.ge2'))
    error('Image types other than .ge2 are not currently supported.')
end

% set default options
if(~isfield(options,'specificDark'))
    options.specificDark = imgInfo.numDark;
end
dome = 4;
debug = false;

% overwrite options if passed in
if(isfield(options,'dome'))
    dome = options.dome(1);
end

if(isfield(options,'debug'))
    debug = options.debug;
end

% find parameters
numLS = size(imgInfo.sweepNumArray,1);
numSweeps = size(imgInfo.sweepNumArray,2);
numGrains = size(allLoadData(1).grainData,2);

% initialize storage variables
pfData = struct();


% for each loadstep
for ii = 1:numLS
    
    grainData = allLoadData(ii).grainData;
    
    % Track the number of peaks used for each grain while processing each
    % loadstep
    numPeaksPerGrain = ones(numGrains,1);
    
    % for each image stack
    for jj = 1:numSweeps
        
        % number of images to skip to load the end of a sweep
        num2skip = imgInfo.numDark + imgInfo.numPerStack - (dome*2);
        
        %+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        %Load the image stack
        if(jj == 1) % the first sweep
            disp(['Loading first sweep in load step #' num2str(ii)])
            
            % for the first fastsweep, load only the first fastsweep and
            % the first 2*dome of the next fastsweep
            
            imgStack = FastsweepRamLoad(imgInfo.dir,imgInfo.stem,imgInfo.sweepNumArray(ii,1),...
                imgInfo.numDark,imgInfo.numPerStack,instr,options,[]);
            
            imgStack = cat(3,imgStack, FastsweepRamLoad(imgInfo.dir,imgInfo.stem,...
                imgInfo.sweepNumArray(ii,2),imgInfo.numDark,dome*2,instr,options,[]) );
            
            
            % set peak search bounds for the sweep
            omeLow = omeBnds(1,1) + (dome+1) * omeStep;
            omeHigh = omeBnds(1,2);
            
            
        elseif(jj == numSweeps) % the last sweep
            disp(['Loading sweep #' num2str(numSweeps) ' in load step #' num2str(ii)])
            
            % for the last fastsweep load the last 2*dome of the previous
            % sweep and the last fastsweep
            
            imgStack = FastsweepRamLoad(imgInfo.dir,imgInfo.stem,imgInfo.sweepNumArray(ii,numSweeps-1),...
                num2skip,dome*2,instr,options,[]);
            
            imgStack = cat(3,imgStack,FastsweepRamLoad(imgInfo.dir,imgInfo.stem,...
                imgInfo.sweepNumArray(ii,numSweeps),imgInfo.numDark,imgInfo.numPerStack,instr,options,[]));
            
            
            % set peak search bounds for the sweep
            omeLow = omeBnds(numSweeps,1);
            omeHigh = omeBnds(numSweeps,2) - (dome+1) * omeStep;
            
            
        else  % all sweeps in the middle
            disp(['Loading sweep #' num2str(jj) ' in load step #' num2str(ii)])
            
            % otherwise, load the end of the previous sweep, the current
            % sweep, and the start of the next sweep
            imgStack = FastsweepRamLoad(imgInfo.dir,imgInfo.stem,imgInfo.sweepNumArray(ii,jj-1),...
                num2skip,dome*2,instr,options,[]);
            
            imgStack = cat(3,imgStack,FastsweepRamLoad(imgInfo.dir,imgInfo.stem,...
                imgInfo.sweepNumArray(ii,jj),imgInfo.numDark,imgInfo.numPerStack,instr,options,[]));
            
            imgStack = cat(3,imgStack, FastsweepRamLoad(imgInfo.dir,imgInfo.stem,...
                imgInfo.sweepNumArray(ii,jj+1),imgInfo.numDark,dome*2,instr,options,[]) );
            
            
            % set peak search bounds for the sweep, one frame in on each side
            omeLow = omeBnds(jj,1)+1*omeStep;
            omeHigh = omeBnds(jj,2)-1*omeStep;        
        end
        
        
        %+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        % itterate over each grain
        for kk = 1:numGrains
            
            disp(['Processing grain #' num2str(kk)])
            
            % find all the peaks inside the omega bounds for the sweep
            %%%
            % This section is used because somehow FindOmeRough(frameno)
            % and ome are out of synch
            extremaOme = [min(min(omeBnds)),max(max(omeBnds))];
            tempOme = FindOmeRough(extremaOme(1),extremaOme(2),omeStep,grainData(kk).fitData.frameno);
            index1 = (tempOme >= omeLow) & (tempOme <= omeHigh);
            %%%
            % This doesn't work for some reason.
            % index1 = (grainData(kk).fitData.ome >= omeLow) & (grainData(kk).fitData.ome <= omeHigh);
            
            if debug
                keyboard
            end
            
            % convert pksToUse to a logical array
            index2 = zeros(size(grainData(kk).fitData.ome,1),1);
            index2(grainData(kk).fitData.pksToUse) = 1;
            index2 = logical(index2);
            
            frameno = grainData(kk).fitData.frameno(index1 & index2);
            pixCoords = grainData(kk).fitData.pixCenter(index1 & index2,:);
            tths = grainData(kk).fitData.tth(index1 & index2);
            hkls = grainData(kk).fitData.hkl(index1 & index2,:);
            
            % shift the global frame number to a frame number in the sweep
            if(jj~=1)
                frameno = frameno - (jj-1)*imgInfo.numPerStack + dome*2;
            end
            
            % integrate each peak
            for ll = 1:size(hkls,1)
                
                % create local storage variable
                peak = struct();
                
                peak.IList = [];
                peak.tthList = [];
                peak.etaList = [];
                peak.omeList = [];
                peak.gS = [];
                peak.hkl = hkls(ll,:);
                
                %%% DEBUG %%%
                debugLogical1 = ismember(peak.hkl,grainData(kk).fitData.hkl(index1 & index2,:));
                if(~debugLogical1)
                    disp([ii, jj, kk, ll])
                    disp(peak.hkl)
                    error('Selected HKL does not exist.')
                end
                %%%%%%%%%%%%%
                
                % integrate each frame
                for mm = frameno(ll)-dome:frameno(ll)+dome
                    
%                     if (mm<=0)
%                         disp(mm)
%                         disp(ll)
%                     else
                    
                    integratedPeak = IntegratePeak(imgStack(:,:,mm),pixCoords(ll,:),instr,cakeParms,options);
                    
                    % build lists to calculate gS with
                    etaList = integratedPeak.etaList;
                    
                    % assign the unstrained tth to all points in the peak
                    tthList = ones(length(etaList),1)*tths(ll);
                    
                    % Compensate for the inclusion of extra frames for
                    % finding omeR
                    if jj==1
                        omeR = FindOmeRough(omeBnds(jj,1),omeBnds(jj,2)+(dome*2)*omeStep,omeStep,mm);
                    elseif jj == numSweeps
                        omeR = FindOmeRough(omeBnds(jj,1)-(dome*2)*omeStep,omeBnds(jj,2),omeStep,mm);
                    else
                        omeR = FindOmeRough(omeBnds(jj,1)-(dome*2)*omeStep,omeBnds(jj,2)+(dome*2)*omeStep,omeStep,mm);
                    end
                    

                    omeList = ones(length(etaList),1)*omeR;
                    
                    qList = qFromTthEtaOme(instr,tthList,integratedPeak.etaList,omeList);
                    
                    % write to storage variable
                    peak.IList = [peak.IList; integratedPeak.etaIList];
                    peak.tthList = [peak.tthList; tthList];
                    peak.etaList = [peak.etaList; etaList];
                    peak.omeList = [peak.omeList; omeList];
                    peak.gS = [peak.gS; qList];
%                     end
                end
                
                % save the data to the output structure
                pfData(ii).grains(kk).peaks(numPeaksPerGrain(kk)) = peak;
                
                % old version (pre 4/2016)
%                 pfData.LoadSteps(ii).grains(kk).peaks(numPeaksPerGrain(kk)) = peak;
                
                % increment numPeaksPerGrain
                numPeaksPerGrain(kk) = numPeaksPerGrain(kk)+1;

            end
        end
    end
end







end