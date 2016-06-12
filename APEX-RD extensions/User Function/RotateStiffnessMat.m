function CmatS = RotateStiffnessMat(CmatC, Rsc)
% RotateStiffnessMat: This function takes the single crystal stiffness
% matrix in the crystal coordinate system and the rotation matrix for the
% crystal's orientation and calculates the stiffness matrix in the sample
% coordinate system.
% 
% USAGE: CmatS = RotateStiffnessMat (CmatC, Rsc)
% 
% AUTHOR: Timothy Long
% 
% INPUTS:
%   CmatC is 6 x 6:
%       The single crystal stiffness tensor represented using the Voigt
%       notation relative to the crystal coordinate system.  This notation
%       matches the notation used in XFormStressTrainVT.
% 
%   Rsc is 3 x 3:
%       The rotation matrix describing the coordinate transformation from
%       the crystal coordinate system (c) to the sample coordinate system
%       (s), Rsc * c = s
% 
% OUTPUTS:
%   CmatS is 6 x 6:
%       The single crystal stiffness tensor represented using the Voigt
%       notation relative to the sample coordinate system.
% 
% NOTES:
%   Started 2015_7_2
%
%   The Voigt convention used for the strain vector that matches this
%   stiffness matrix is:
%       strainTen_23 = 1/2 * strainVec_4 
%       strainTen_13 = 1/2 * strainVec_5 
%       strainTen_12 = 1/2 * strainVec_6
%   This matches the convention used in XFormStressStrainVT.  Make sure to
%   use the appropriate elastic constants!


K1 = Rsc.^2;

K2 = [Rsc(1,2)*Rsc(1,3), Rsc(1,3)*Rsc(1,1), Rsc(1,1)*Rsc(1,2);
    Rsc(2,2)*Rsc(2,3), Rsc(2,3)*Rsc(2,1), Rsc(2,1)*Rsc(2,2);
    Rsc(3,2)*Rsc(3,3), Rsc(3,3)*Rsc(3,1), Rsc(3,1)*Rsc(3,2)];

K3 = [Rsc(2,1)*Rsc(3,1), Rsc(2,2)*Rsc(3,2), Rsc(2,3)*Rsc(3,3);
    Rsc(3,1)*Rsc(1,1), Rsc(3,2)*Rsc(1,2), Rsc(3,3)*Rsc(1,3);
    Rsc(1,1)*Rsc(2,1), Rsc(1,2)*Rsc(2,2), Rsc(1,3)*Rsc(2,3)];

K4 = [Rsc(2,2)*Rsc(3,3)+Rsc(2,3)*Rsc(3,2), Rsc(2,3)*Rsc(3,1)+Rsc(2,1)*Rsc(3,3), Rsc(2,1)*Rsc(3,2)+Rsc(2,2)*Rsc(3,1);
    Rsc(3,2)*Rsc(1,3)+Rsc(3,3)*Rsc(1,2), Rsc(3,3)*Rsc(1,1)+Rsc(3,1)*Rsc(1,3), Rsc(3,1)*Rsc(1,2)+Rsc(3,2)*Rsc(1,1);
    Rsc(1,2)*Rsc(2,3)+Rsc(1,3)*Rsc(2,2), Rsc(1,3)*Rsc(2,1)+Rsc(1,1)*Rsc(2,3), Rsc(1,1)*Rsc(2,2)+Rsc(1,2)*Rsc(2,1)];

K = [K1 2*K2; K3 K4];

CmatS = K * CmatC * K';

end