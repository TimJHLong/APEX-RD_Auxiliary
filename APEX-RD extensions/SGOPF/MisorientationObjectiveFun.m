function delta = MisorientationObjectiveFun(w, gS, traj, peakNum, options)
% An objective function to minimize the residual between the predicted
% (from w) trajectories and the measured trajectories as a function of w.

% defaults
checkNeg = false;
zetaWeight = 5e7;
minTrajMag = 1e6; % ~1e4 smaller than gS
maxAngle = 45; % largest permitted angle between measured and calcualted
normMode = 0;

if exist('options','var')
    if isfield(options,'checkNeg')
        checkNeg = options.checkNeg;
    end
    
    if isfield(options,'zetaWeight')
        zetaWeight = options.zetaWeight;
    end
    
    if isfield(options,'minTrajMag')
        minTrajMag = options.minTrajMag;
    end
    
    if isfield(options,'maxAngle')
        maxAngle = options.maxAngle;
    end
    
    if isfield(options,'normMode')
        normMode = options.normMode;
    end
end

numPeaks = size(peakNum,1);

if isinf(norm(w))
    delta = inf;
    return
end

% calculate trajectory from wS for each gS
magG = sum(gS.^2,2).^0.5;
nS = [gS(:,1)./magG, gS(:,2)./magG, gS(:,3)./magG];
vS = [w(2)*nS(:,3)-w(3)*nS(:,2), w(3)*nS(:,1)-w(1)*nS(:,3),...
    w(1)*nS(:,2)-w(2)*nS(:,1)];

% removed vS's that are too small.  This incures no penalty.
magV = sum(vS.^2,2).^0.5;
index1 = magV < minTrajMag;
tempMagV = magV(~index1);
tempVS = vS(~index1,:);
tempTraj = traj(~index1,:);
tempPeakNum = peakNum(~index1);

% check both +traj and -traj
magT = sum(tempTraj.^2,2).^0.5;
if(normMode ~= 1)
    % difference of magnitudes
    u = magT.^2 - tempMagV.^2;
else
    % magnitude of differences
    % not yet working
    deltaT = tempVS - tempTraj;
    u = sum(deltaT.^2,2).^0.5;
end

% check angle
dots = tempTraj(:,1).*tempVS(:,1) + tempTraj(:,2).*tempVS(:,2) + ...
    tempTraj(:,3).*tempVS(:,3);
zeta = acosd(dots./(magT.*tempMagV));
% filter out NaN's
zeta(isnan(zeta))=0;
if checkNeg % if the angle is more than 90, use the conjugate angle
    zeta(zeta>90) = 180 - zeta(zeta>90);
end

% removed trajectories with too large of angles
index2 = zeta > maxAngle;
zeta(index2) = [];
u(index2) = [];
tempPeakNum2 = tempPeakNum(index2);
zeta = zeta * zetaWeight;

% check that all peaks with measurable trajectories are used
resid = 0;
for ii = 1:numPeaks
    check1 = ismember(ii,tempPeakNum2);
    [check2,loc] = ismember(ii,tempPeakNum);
    % if none of the trajectories in a peak match the calculated, incure a
    % penalty equal to the average magnitude of the trajectories in that
    % peak.
    if ~check1 && check2
        tempT = tempTraj(loc,:);
        tempMag = sum(tempT.^2,2).^0.5;
        resid = resid + sum(tempMag)/length(tempMag);
    end
end
delta = sum(u.^2 + zeta.^2) + resid^2;
if isnan(delta)
    keyboard
end
end