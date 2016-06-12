function [strainTenS, strainVecS, finalLatticeStrain, finalNormalS,...
    index] = RefineStrainTensor(fitData0, fitData, options)
% RefineStrainTensor
%   This function calculates the strain tensor from change in lattice
%   strains and has several options to filter the lattice strains for a
%   better fit
% 
% USAGE: 
%   [strainTenS, strainVecS, finalLatticeStrain, finalNormalS] = ...
%   RefineStrainTensor(fitData0, fitData, options);
% 
% AUTHOR: Timothy Long 
% 
% INPUTS:
%   fitData0 is a 1x1 structure
% 
%   fitData is a 1x1 structure
% 
%   options is a 1x1 structure
% 
% OUTPUTS:
%   strainTenS is 3x3:
%       The strain tensor in the sample coordinate system represented as a
%       3x3 matrix.
% 
%   strainVecS is 6x1
%       The strain tensor in the sample coordinate system represented in
%       the Voigt notation.
% 
%   finalLatticeStrain is n (number of planes) x 1
%       The lattice strains that were used to calculate the strain tensor
%
%   finalNormalS is n x 1
%       The normalized gS0 vectors corresponding to the lattice planes that
%       were used to calculate the strain tensor
%
%   index is n x 1
%       The indecies along the first dimension of the gS0 (after filtering)
%       of the peaks used in the final strain tensor.  This is NOT Miller
%       Indicies.
%
% ACKNOWLEDGEMENTS:
%   Darren Pagan and Mark Obstalecki (T.L. uses their code as references)
%
% NOTES:
%   Started 2015_Aug_28
%
%


% set default options
pksToUse=1:min(length(fitData0.tth),length(fitData.tth));
filterStrains = false;
matchPks = false;
useGs = false;
itterate = true;
latticeStrainTol = 1e-4;
deltaStrainTol = 1e-4;
numLoops = 10;
debug = false;


% overwrite options if passed in
if(exist('options','var'))
    
    if(isfield(options,'RwThreshold'))
        RwThreshold = options.RwThreshold;
        pksToUse=GoodFitPeaks(fitData,RwThreshold);
        pksToUse=pksToUse(pksToUse<=length(fitData0.tth));
    end
    
    if(isfield(options,'pksToUse'))
        pksToUse = options.pksToUse;
    end
    
    if(isfield(options,'filterStrains') && isfield(options,'maxStrains'))
        filterStrains = options.filterStrains;
        maxStrains = options.maxStrains;
    end
    
    if(isfield(options,'matchPks'))
        matchPks = options.matchPks;
    end
    
    if(isfield(options,'useGs'))
        useGs = options.useGs;
    end
    
    if(isfield(options,'itterate'))
        itterate = options.itterate;
    end
    
    if(isfield(options,'latticeStrainTol'))
        latticeStrainTol = options.latticeStrainTol;
    end
    
    if(isfield(options,'deltaStrainTol'))
        deltaStrainTol = options.deltaStrainTol;
    end
    
    if(isfield(options,'numLoops'))
        numLoops = options.numLoops;
    end
    
    if(isfield(options,'debug'))
        debug = options.debug;
    end
    
end


% calculate lattice strains
if matchPks
    [peakList0, peakList] = FindMatchingPeaks(fitData0, fitData);
    
    tth0 = fitData0.tth(peakList0);
    tth = fitData.tth(peakList);
    gS0 = fitData0.gS(peakList0,:);
    gS = fitData.gS(peakList,:);
else
    tth0 = fitData0.tth(pksToUse);
    tth = fitData.tth(pksToUse);
    gS0 = fitData0.gS(pksToUse,:);
    gS = fitData.gS(pksToUse,:);
end
nS0 = UnitVector(gS0')';

if useGs
    % This is used for correcting energy from stress
    latticeStrain = (sum(gS0.^2,2).^0.5)./(sum(gS.^2,2).^0.5)-1;
else
    % This is the standard method
    latticeStrain = sind(tth0/2)./sind(tth/2) - 1;
end

if(filterStrains)
%   filter out strains that are larger than maxStrains(1) and smaller than
%   maxStrains(2)
    lIndex1 = latticeStrain > maxStrains(1) | latticeStrain < maxStrains(2);
    
    latticeStrain(lIndex1)=[];
    nS0(lIndex1,:) = [];  
end

%Calculate Strain Tensor components in the sample frame
strainVecS = StrainRosette(latticeStrain,nS0);
strainTenS= XFormStressStrainVT(strainVecS,'strain');

% record the lattice strains used and the associated lattice plane normals
finalLatticeStrain = latticeStrain;
finalNormalS = nS0;


if(itterate)
    % initialize loop variables
    continueLoop = true;
    loopTally = 1;
    
    while continueLoop
        % calculate the lattice strains from the fit strain tensor
        calcLS = zeros(size(nS0,1),1);
        
        for ii=1:size(nS0,1)
            calcLS(ii) = nS0(ii,:) * strainTenS * nS0(ii,:)';
        end
        
        lsError = abs(calcLS - latticeStrain);
        
        % filter out lattice strains that are too far from the fit
        index = lsError < latticeStrainTol;
        tempLS = latticeStrain(index);
        tempNS = nS0(index,:);
        
        % calculate new strain tensor
        tempStrainVS = StrainRosette(tempLS,tempNS);
        tempStrainTenS = XFormStressStrainVT(tempStrainVS,'strain');
        
        % check to see if the strain tensor changed much
        deltaStrainTS = tempStrainTenS - strainTenS;
        magDStrainTS = sqrt( trace(deltaStrainTS^2) );
        
        % if the change was less than the tolerance, exit the loop
        if magDStrainTS < deltaStrainTol || loopTally >= numLoops
            continueLoop = false;
        else
            loopTally = loopTally+1;
        end
        
        % write tempoarary loop variables to permanent ones
        strainTenS = tempStrainTenS;
        strainVecS = tempStrainVS;  
        finalLatticeStrain = tempLS;
        finalNormalS = tempNS;
    end  
end

if(debug)
    disp(['Number of Itterations: ' num2str(loopTally)])
    disp(['Magnitude of delta: ' num2str(magDStrainTS)])
    disp(['Number of peaks: ' num2str(length(finalLatticeStrain))])
    disp('-----')
    disp(' ')
end

end