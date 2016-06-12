function [finalWS,finalResids,finalFlags] = MultiTrajCalcMisorientationAxis(peakData, options)
% MultiTrajCalcMisorientationAxis: This function takes the trajectory data
% from all peaks for one grain in one load step and attempts to fit a
% misorientation vector to one trajectory from each peak.  Optionally,
% additional misorientation vectors can be calculated, up to the number of
% trajectories per peak.
% 
% USAGE: 
%   [finalWS,finalResids,finalFlags] =...
%       MultiTrajCalcMisorientationAxis(trajData,options);
%
% AUTHOR: Timothy Long
% 
% CITATION:
%   Pagan, D.  "Determining Heterogeneous Slip Activity on Multiple Slip
%       Systems from Single Crystal Orientation Pole Figures."
% 
% INPUTS:
%   peakData is a 1 x n structure:
%       A structure containing the trajectory data for each peak from one
%       grain.  Fields are:
%       
%       .trajStart is 3 x n (# of trajectories):
%           A list of column vectors giving the starting points of all n
%           trajectories for one peak.  Units are m^-1.
%       .trajEnd is 3 x n:
%           A list of column vectors giving the ending points of all n
%           trajectories for one peak.  Units are m^-1.
%       .trajDir is 3 x n:
%           A list of unit column vectors giving the direction of all n
%           trajectories for one peak.  Unitless.
%       .trajMag is n x 1:
%           A list of the magnitudes of each trajectory.  Units are m^-1.
%       .gS is 1 x 3:
%           The reciprocal lattice vector of the hkl planes that generated
%           the peak, represented in the sample coordinate system.  Units
%           are m^-1.
% 
%   options is a 1 x 1 structure;
%       A structure containing options than change how the code executes.
%       Fields used are:
%       
%       .numAxes is 1 x 1:
%           The number of misorientation axes to try and find.  Maximum is
%           n, the number of trajectories per peak.  Default is 1.
%       .peakFitFilter is a boolean
%           If true, filters the peaks by the Rw values from peak fitting.
%           If false, uses all peaks.  Default is false.
%       .RwThreshold is 1 x 3
%           The Rw values used to filter out poorly fit peaks.  Peaks with
%           an Rw value greater than one of the thresholds are removed from
%           the peak & trajectory lists.  Order is RwE,RwT,RwO.  Only
%           needed if options.peakFitFilter = true
%       .maxFEval is 1 x 1
%           The maximum number of function evalutations the minimization is
%           allowed to make.  Default is 2000.
%       .maxIt is 1 x 1
%           The maximum number of itterations the minimization is allowed
%           to make.  Default is 2000.
% 
% OUTPUTS:
%   wS is 3 x m (# of misorientation axes):
%       A list of column vectors giving the calculated misorientation axes
%       in the sample coordinate system.
%
%   fitQ is 1 x m:
%       A metric of the quality of the fit of each calculated trajectory to
%       the data.

% defaults
numAxes = 1;
peakFitFilter = false;
maxFEval = 2000;
maxIt = 2000;
debug = false;

% overwrite defaults if passed in
if exist('options','var')
    if isfield(options,'numAxes')
        numAxes = options.numAxes;
    end
    
    if isfield(options,'peakFitFilter')
        peakFitFilter = options.peakFitFilter;
    end
    
    if isfield(options,'RwThreshold')
        RwThreshold = options.RwThreshold;
    else
        % peak filtering requires having RwThresholds
        peakFitFilter = false;
    end
    
    if isfield(options,'maxFEval')
        maxFEval = options.maxFEval;
    end
    
    if isfield(options,'maxIt')
        maxIt = options.maxIt;
    end
    
    if isfield(options,'debug')
        debug = options.debug;
    end
end

numTs = size(peakData(1).trajDir,2);
numPeaks = size(peakData,2);

if numAxes > numTs
    numAxes = numTs;
end

% build input variables
gS = [peakData.gS]';
traj = [peakData.traj]';
peakNum = repmat((1:numPeaks),numTs,1);
peakNum = peakNum(:);

% filter out peaks with low Rw values
if peakFitFilter
    % reformat data
    disp('Filtering peaks with poor fits')
    tempFitData = struct();
    tempFitData.RwE = [peakData.RwE];
    tempFitData.RwR = [peakData.RwR];
    tempFitData.RwO = [peakData.RwO];
    goodPeaks = GoodFitPeaks(tempFitData,RwThreshold);
    
    % apply filter
    gS = gS(goodPeaks,:);
    traj = traj(goodPeaks,:);
    peakNum = peakNum(goodPeaks,:);
end

% pre-allocate storage variables
wS = zeros(3,numTs);
resids = zeros(1,numTs);
flags = zeros(1,numTs);

% guess the misorientation axes
searchBnds = norm(traj(1,:));
wSGuesses = RoughFindWS(traj', gS', numTs, searchBnds);

% create objective function
solverOpts = optimset('MaxFunEvals',maxFEval,'MaxIter',maxIt);
fun = @(wS) MisorientationObjectiveFun(wS, gS, traj, peakNum, options);
for ii = 1:size(wSGuesses,2)
    % guess the misorientation axis
    % build guess from one trajectory
%     v1 = peakData(1).trajStart(:,ii);
%     v2 = peakData(1).trajEnd(:,ii);
%     wGuess = cross(v1,v2)/norm(v1);

    wGuess = wSGuesses(:,ii)';
    % I should probably make everything consistent so I don't need all
    % these transposes floating around.
    
    % fit misorientation axis to all trajectories
    [wS(:,ii),resids(ii),flags(ii)] = fminsearch(fun,wGuess,solverOpts);
    
    % remove used trajectories (not yet implemented)
end

[resids,index1] = sort(resids);
wS = wS(:,index1);
flags = flags(index1);

% return best axes
finalWS = wS(:,1:numAxes);
finalResids = resids(1:numAxes);
finalFlags = flags(1:numAxes);

if debug
    disp(['Residuals: ' num2str(resids')])
    disp(['Exit Flags:' num2str(flags')])
    disp('Misorientation axes:')
    disp(num2str(wS))
    keyboard
end
end
