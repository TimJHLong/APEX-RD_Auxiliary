function trajData = MainSGPFAnalysis(pfData, allLoadData, threshold, grainNum, instr, options)
% MainSGPFAnalysis: This function takes the pole figure data of one grain
% from any number of load steps and finds the trajectories of all (or a
% user defined list) hkl's at all loads.  A misorientation vector is fit to
% all the trajectories.  Optionally, the misorientation vector is
% decomposed into relative slip activity.
% 
% USAGE: 
%   [trajData, slipData] = MainSGPFAnalysis(pfData, allLoadData,...
%       threshold, omeStep, instr, material, options);
% 
% AUTHOR: Timothy Long
%
% INPUTS:
%   pfData is a 1 x ls (# of load steps) structure
%       fields are:
%       .grains is a 1 x ng (# of grains) structure
%           fields are:
%           .peaks is a 1 x np (# of integrable peaks) structure
%           fields are:
%               .omeList is n x 1:
%                   A list of all the omega coordinates of all the points
%                   in the pole figure.
%               .etaList is n x 1:
%                   A list of all the eta coordinates of all the points in
%                   the pole figure.
%               .IList is n x 1:
%                   A list of the two-theta integrated intensity of each
%                   point in the pole figure.
% 
%   allLoadData is a 1 x ls structure
%       fields are:
%       .grainData is a 1 x ng structure
%           needed fields are:
%           .fitData is a 1 x 1 structure
%           needed fields are:
%               .ome is m (# of peaks) x 1:
%                   The omega position of all m peaks
%               .tth is m x 1:
%                   The two theta position of all m peaks
%               .eta is m x 1:
%                   The eta position of all m peaks
%
%   threshold is a scalar:
%       A threshold used to convert the pole figure to a binary image for
%       processing.  Points with intensities below the threshold are set to
%       0, and points equal to or above it are set to 1.  The value is a
%       percentage relative to the intensity of the maxima closest to the
%       fit peak center.
% 
%   options is a structure
%       The options structure.  Used fields are:
%
% OUTPUTS:
%   trajData is a 1 x ls structure:
%       Contains data describing the trajectories and misorientation axes
%       for one grain in all loadsteps. Fields are:
%       
%       .peaks is a 1 x m structure:
%           Contains the trajectory information for each peak.  Has fields:
%       
%           .trajStart
% 
%       .wS is 3 x k (# of misorientation axes)
%           A list of k column vectors giving the best fit misorientation
%           axes, sorted by lowest residual
% 
%       .residuals is k x 1:
%           The residuals for each misorientation axis.
%
%       .flags is k x 1:
%           The exit flags from the solver for each wS


%
% Restructure pf data
%
tempPFData = struct();

for ii = 1:size(allLoadData,2)
    tempPFData(ii).peaks = pfData(ii).grains(1).peaks;
end

%
% Normalize pole figures by max intensity
%

for ii = 1:size(tempPFData,2)
    for jj = 1:size(tempPFData(ii).peaks,2)
        tempPFData(ii).peaks(jj).IList = 100*tempPFData(ii).peaks(jj).IList/max(tempPFData(ii).peaks(jj).IList);
    end
end

%
% Automatic calculation of the pf trajectory
%

trajData = struct();
for ii = 1:size(tempPFData,2)
    disp(['Processing Grain #' num2str(grainNum) ' in load step #' num2str(ii)])
    
    hkls = allLoadData(ii).grainData(grainNum).fitData.hkl;
    
    for jj = 1:size(tempPFData(ii).peaks,2)
        % find the peak in allLoadData
        [~,index] = ismember(tempPFData(ii).peaks(jj).hkl,hkls,'rows');
        
        % make sure fitEta and the pole figure integration window use the
        % same eta range
        tempPFData(ii).peaks(jj).eta0 =...
            allLoadData(ii).grainData(grainNum).fitData.eta(index);
        if tempPFData(ii).peaks(jj).eta0 > 180
            tempPFData(ii).peaks(jj).eta0 = tempPFData(ii).peaks(jj).eta0 - 360;
        elseif tempPFData(ii).peaks(jj).eta0 < -180
            tempPFData(ii).peaks(jj).eta0 = tempPFData(ii).peaks(jj).eta0 + 360;
        end
        
        % write fit omega, two theta, and gS values to the temporary data
        % structure.
        tempPFData(ii).peaks(jj).ome0 =...
            allLoadData(ii).grainData(grainNum).fitData.ome(index);
        tempPFData(ii).peaks(jj).tth0 =...
            allLoadData(ii).grainData(grainNum).fitData.tth(index);
        tempPFData(ii).peaks(jj).gS0 = qFromTthEtaOme(instr,...
            tempPFData(ii).peaks(jj).tth0, tempPFData(ii).peaks(jj).eta0,...
            tempPFData(ii).peaks(jj).ome0)';
        
        % calculate the trajectory of the peak
        [tempStart, tempEnd, trajData(ii).peaks(jj).trajMag] =...
            AutoCalcPFTrajectory(tempPFData(ii).peaks(jj),...
            threshold, tempPFData(ii).peaks(jj).eta0,...
            tempPFData(ii).peaks(jj).ome0, 1/20, options);
        
        % Convert from (eta,ome) to gS
        tth = repmat(tempPFData(ii).peaks(jj).tthList(1),options.numTrajs,1);
        trajStart = qFromTthEtaOme(instr,tth,tempStart(1,:)',tempStart(2,:)');
        startMag = sum(trajStart.^2,2).^0.5;
        trajEnd = qFromTthEtaOme(instr,tth,tempEnd(1,:)',tempEnd(2,:)');
        endMag = sum(trajEnd.^2,2).^0.5;
        
        % calculated values
        traj = trajEnd - trajStart;
        trajMag = sum(traj.^2,2).^0.5;
        trajDir = traj./repmat(trajMag,1,3);
        
        % store data
        trajData(ii).peaks(jj).trajStart = trajStart';
        trajData(ii).peaks(jj).trajSMag = startMag;
        trajData(ii).peaks(jj).trajEnd = trajEnd';
        trajData(ii).peaks(jj).trajEMag = endMag;
        trajData(ii).peaks(jj).trajDir = trajDir';
        trajData(ii).peaks(jj).trajMag = trajMag;
        trajData(ii).peaks(jj).traj = traj';
        trajData(ii).peaks(jj).gS =...
            repmat(tempPFData(ii).peaks(jj).gS0,1,size(trajStart,1));
        
        % copy peak fit quality to trajData for ease of use
        trajData(ii).peaks(jj).RwE =...
            allLoadData(ii).grainData(grainNum).fitData.RwE(index);
        trajData(ii).peaks(jj).RwR =...
            allLoadData(ii).grainData(grainNum).fitData.RwR(index);
        trajData(ii).peaks(jj).RwO =...
            allLoadData(ii).grainData(grainNum).fitData.RwO(index);
    end
   % find misorientation axes
   [wS,residuals,flags] =...
       MultiTrajCalcMisorientationAxis([trajData(ii).peaks], options);
   trajData(ii).wS  = wS;
   trajData(ii).residuals = residuals;
   trajData(ii).flags = flags;
end

end