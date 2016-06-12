function qS = PrecessionCorrectQ(tth, eta ,ome ,pS, instr)
% PrecessionCorrectQ:  Takes measured tth, eta, and ome values of
% diffraction events for one grain, plus the best fit position vector (in
% the sample coorinates) of that grain and calculates the true scattering
% vectors by correcting for precession.
% 
% USAGE: qS = PrecessionCorrectQ(tth, eta ,ome ,pS, instr);
%
% AUTHOR: Timothy Long 
% 
% INPUTS:
%   tth is n (# of peaks) x 1:
%       A list of the measured two theta values for each diffraction event.
% 
%   eta is n x 1:
%       A list of the measured eta values for each diffraction event.
% 
%   ome is n x 1:
%       A list of the omega values for each diffraction event.
%
%   pS is 3 x 1:
%       The position of the center of mass of the grain in the sample
%       coordinate system.
%
%   instr is a structure:
%       The instrument data structure.
%
% OUTPUT:
%   qS is n x 3:
%       The precession corrected scattering vectors of each input peak.
%
% NOTES:
%   Use "PrecessionFitObjective.m" as a reference.

% calculate measured position from input variables
d = repmat(instr.distance,size(tth));
r_mes = d .* tand(tth);
X_mes = r_mes .* cosd(eta);
Y_mes = r_mes .* sind(eta);

% calculate pL at each diffraction event
pL = zeros(3,size(tth,1));
for ii = 1:size(tth,1)
    Rls=[cosd(ome(ii)),0,sind(ome(ii));0,1,0;-sind(ome(ii)),0,cosd(ome(ii))];
    pL(:,ii) = Rls * pS;
end

% adjust for position
d_corr = d + pL(3,:)'; % true distance
X_corr = X_mes - pL(1,:)'; % peak X relative to grain center
Y_corr = Y_mes - pL(2,:)'; % peak Y relative to grain center
tth_corr = 180/pi * atan((X_corr.^2 + Y_corr.^2).^0.5 ./ d_corr);
eta_corr = 180/pi * atan2(Y_corr,X_corr);

% calculate position subtracted corrected gS
qS = qFromTthEtaOme(instr,tth_corr,eta_corr,ome);

%%% DEBUG %%%
% keyboard
%%%%%%%%%%%%%
end