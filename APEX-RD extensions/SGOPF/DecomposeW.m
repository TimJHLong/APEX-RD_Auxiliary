function slipActivity = DecomposeW(w, stressTenC, material, options)
% DecomposeW: This function takes a misorientation vector, w, the stress
% tensor in the crystal frame, and the slip systems of the crystal and
% attempts to break w down into slip activity on the slip systems with the
% largest resolved shear.  This function is only valid in the limit of
% small plastic strains and a final grain rotation of approximately 0.
% 
% USAGE: 
%   slipActivity = DecomposeW(w, stressTenC, material, options);
% 
% AUTHOR: Timothy Long
% 
% CITATION:
%   Pagan, D. C. Miller, M. P.  "Determining Heterogeneous Slip Activity on
%   Multiple Slip Systems from Single Crystal Orientation Pole Figures."
%   Publication Pending.
% 
% INPUTS:
%   w is 1 x 3
%       The misorientation vector of the grain
% 
%   stressTenC is 3 x 3
%       The stress tensor of the grain, represented in matrix form in the
%       crystal coordinate system.
% 
%   material is a structure
%       The material structure.  Used fields are:
%       
%       .slip is a structure
%           Contains fields that describe the slip systems of the crystal
% 
%           .burgersVecs is n x 3:
%               The burgers vectors for all 'n' possible slip systems.
%               May contain repeated vectors.
%
%           .planeNorms is n x 3:
%               The slip plane normals corresponding with the burgers
%               vectors in material.slip.burgersVecs.  
%
%   options is a structure
%       The options structure
% 
% OUTPUTS:
%   slipActivity is n x 1:
%       Relative slip system activity as calculated by decomposing w.
%       slipActivity is normalized to sum to 1.
% 
% NOTES:
%   Started 2016/2/3
%
%   The calculations in this function are based on the citated paper and
%   make several major approximations:
%       1: A viscoplastic model is used to approximate F
%       2: The total rotation of the grain, R, is assumed to be negligable,
%          i.e. R nearly equals I (eye(3))
%       3: The plastic strains are small and thus can be linearly
%          decomposed.
%       4: The relative activity of each slip system is constant over the
%          grain.  I.e. the pole figures are smeared along arcs of a circle
%          on the sphere, not some other shape.

% default options
debug = false;

% overwrite defaults
if exist('options','var')
    if isfield(options,'debug')
        debug = options.debug;
    end
end

% read in slip systems from material structure
s = material.slip.burgersVecs;
n = material.slip.planeNorms;
numSys = size(s,2);

% ensure s & n are unit vectors
s = UnitVector(s')';
n = UnitVector(n')';

% Calculate the rotation axis (m), Schmid Tensor, and resolved shear for
% each slip system.
m = zeros(numSys,3);
schmidTenC = zeros(3,3,numSys);
tau_resolved = zeros(numSys,1);
for ii=1:numSys
    schmidTenC(:,:,ii) = s(ii,:)' * n(ii,:);
    tau_resolved(ii) = abs(trace(schmidTenC(:,:,ii) * stressTenC));
    m(ii,:) = cross(s(ii,:),n(ii,:));
end

% sort slip systems by magnitude of resolved shear
[tauSorted,index1] = sort(tau_resolved);
mSorted = m(index1,:);

% Use the five highest resovled shear slip systems.  If there are not five
% systems with tau~=0, then only use the non-zero ones.  Five slip systems
% are needed for general deformation.
mTemp = mSorted(1:5,:);
mTemp(tauSorted(1:5)==0) = [];
x0 = [norm(w);zeros(size(mTemp,1)-1,1)];

fun = @(x) ObjectiveFunction(x, mTemp,w);
[deltaGamma,check,flag] = fminsearch(fun,x0);

slipActivity = deltaGamma./sum(deltaGamma);

if debug
    disp('Relative slip activity:')
    disp(num2str(slipActivity))
    disp(['Residual = ' num2str(check)])
    disp(['Exit flag = ' num2str(flag)])
    disp('----------')
    disp(' ')
end

end


function check = ObjectiveFunction(deltGamma, m, w)
% A helper function that takes a guess at deltaGamma, a list of the m
% vectors for each slip system, and the misorientation vector, w, and
% calculates the norm of the difference between the measured w and the w
% calculated from deltaGamma and the m's.

test = sum(deltGamma/2.*m,1);
check = norm(w-test);
end