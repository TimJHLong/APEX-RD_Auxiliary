function wS = CalcMisorientationAxis(gS, trajDirs, trajMag, options)
% CalcMisorientationAxis: Using a list of reciprocal lattice vectors and
% the pole figure trajectories assosiated with each gS and calculates the
% best fit misorientation vector, wS ( = delta_phi * tS)
% 
% USAGE: wS = CalcMisorientationAxis(gS, trajDir, trajMag, options);
% 
% AUTHOR: Timothy Long
% 
% SOURCE: 
%   Pagan, D.  "Determining Heterogeneous Slip Activity on Multiple Slip
%   Systems from Single Crystal Orientation Pole Figures."
% 
% INPUTS:
%   gS is n x 3:
%       A list of the reciprocal lattice vectors for one grain
% 
%   trajDir is n x 3:
%       A list of the trajectory directions of the orientation pole
%       assosiated with each gS.
% 
%   trajMag is n x 1:
%       A list of the magnitude of the trajectories of the orientation pole
%       assosiated with each gS.
%
%   options is a structure
%       The options structure
% 
% OUTPUT:
%   w is 1 x 3:
%       The best fit misorientation vector, delta_phi * tS
% 
% NOTES:
%   Started 2016/2/1

% set default options
debug = false;

if exist('options','var')
    if isfield(options,'debug')
        debug = options.debug;
    end
end

% calculate normal vectors
magG = sum(gS.^2,2);
nS = gS./[magG, magG, magG];
traj = trajDirs.*repmat(trajMag,1,3);

fun = @(w) trajectoryObjectiveFun(w, nS, traj);

% [wS, check, ~, flag] = lsqnonlin(fun,([3.5659, 7.1319, 3.5659]*1e8));
[wS, check, flag] = fminsearch(fun,([3.5659, 7.1319, 3.5659]*1e8));

if debug
    disp('Best fit misorientation vector:')
    disp(num2str(wS))
    disp(['Residual = ' num2str(check)])
    disp(['Exit flag: ' num2str(flag)])
end

end


function deltas = trajectoryObjectiveFun(w, nS, traj)
% An objective function to minimize the residual between the predicted
% (from w) trajectories and the measured trajectories as a function of w.
vS = [w(2)*nS(:,3)-w(3)*nS(:,2), w(3)*nS(:,1)-w(1)*nS(:,3),...
    w(1)*nS(:,2)-w(2)*nS(:,1)];

% check both +traj and -traj
% use difference of the norms instead of norm of the difference because I
% also check the angle.
magV = sum(vS.^2,2).^0.5;
magT = sum(traj.^2,2).^0.5;
u = magT.^2 - magT.^2;
dots = traj(:,1).*vS(:,1) + traj(:,2).*vS(:,2) + traj(:,3).*vS(:,3);
zeta = acosd(dots./(magT.*magV));
zeta(zeta>90) = 180 - zeta(zeta>90);
zeta = zeta * 5e7; % weighting factor so both parts have similar magnitude.

% deltas = u.^2 + zeta.^2; % for lsqnonling
deltas = sum(u.^2 + zeta.^2); % for fminsearch
end








