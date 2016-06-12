function [wSGuesses, accNum] = RoughFindWS(trajs, gS, numTrajs, searchBnds)
% RoughFindWS:  This function looks for accumulation points in wS space to
% generate good starting points for the minimizer.
% 
% USAGE: wSGuesses = RoughFindWS(trajs, gS, numTrajs, searchBnds);
%
% AUTHOR: Timothy Long
%
% INPUTS:
%   trajs is 3 x n (total # of trajectories):
%       A list of n trajectories.  Needs at least 2 trajectories to work
% 
%   gS is 3 x n:
%       The reciprocal lattice vectors for each peak
%
%   numTrajs is a scalar
%       The number of trajectories per peak
%
%   searchBnds is 3 x 1
%       The distance away from a point to look for other points
% 
% OUTPUTS:
%   wSGuesses is 3 x m
%       A list of the accumulation points, sorted by number of lines in the
%       accumulation point.  'm' is at least one and at most 'numTrajs'.
%
%   accNum is m x 1:
%       A list of the number of points within +-searchBnds of each of the
%       wSGuesses
% 
% NOTES:
%   Started 2016/Mar/16
%
%   This script takes advantage of the equation vS = wS x nS, where vS is
%   the trajectory of a peak, wS is a misorientation vector of a grain, and
%   nS is the unit vector in the direction of the reciprocal lattice vector
%   of the peak.  Therefore, all trajectories created by one misorientation
%   vector are all orthogonal to that misorientation vector.
%
%   All points will have a numPts >= 1 because they always count
%   themselves.

if length(searchBnds)==1
    searchBnds = searchBnds*ones(3,1);
end

gMag = sum(gS.^2,1).^0.5;
nS = [gS(1,:)./gMag; gS(2,:)./gMag; gS(3,:)./gMag];

% generate all possible wS vectors
wSList = [];
for ii = 1:size(trajs,2)
    a = trajs(:,ii);
    n = nS(:,ii);
    
    % start on the next peak
    temp = mod(ii,numTrajs);
    if temp==0
        temp = numTrajs;
    end
    startInd = ii + (numTrajs-temp+1);
    
    b = trajs(:,startInd:end);
    m = nS(:,startInd:end);
    d = [n(1)./m(1,:); n(2)./m(2,:); n(3)./m(3,:)];

    w1 = (a(3) - b(3,:).*d(1,:))./(n(2) - m(2,:).*d(1,:));
    w2 = (a(1) - b(1,:).*d(2,:))./(n(3) - m(3,:).*d(2,:));
    w3 = (a(2) - b(2,:).*d(3,:))./(n(1) - m(1,:).*d(3,:));
    wSList = cat(2,wSList,[w1;w2;w3]);
    
%     for jj = startInd:size(trajs,2)
%         b = trajs(:,jj);
%         m = nS(:,jj);
%         d = n./m;
%         
%         w1 = (a(3) - b(3)*d(1))/(n(2) - m(2)*d(1));
%         w2 = (a(1) - b(1)*d(2))/(n(3) - m(3)*d(2));
%         w3 = (a(2) - b(2)*d(3))/(n(1) - m(1)*d(3));
%         wSList = cat(2,wSList,[w1;w2;w3]);
%     end
end

% find clusterings of accumulation points
numPts = zeros(size(wSList,2),1);
for ii = 1:size(wSList,2)
    w = wSList(:,ii);
    check1P = wSList(1,:) <= w(1)+searchBnds(1);
    check1N = wSList(1,:) >= w(1)-searchBnds(1);
    check1 = check1P & check1N;
    
    check2P = wSList(2,:) <= w(2)+searchBnds(2);
    check2N = wSList(2,:) >= w(2)-searchBnds(2);
    check2 = check2P & check2N;
    
    check3P = wSList(3,:) <= w(3)+searchBnds(3);
    check3N = wSList(3,:) >= w(3)-searchBnds(3);
    check3 = check3P & check3N;
    
    check = check1 & check2 & check3;
    numPts(ii) = sum(check);
end
[sortedNumP,index] = sort(numPts,1,'descend');
sortedWS = wSList(:,index);

% call center point of largest clusters wSGuesses
wSGuesses = zeros(3,numTrajs);
accNum = zeros(numTrajs,1);

wSGuesses(:,1) = sortedWS(:,1);
accNum(1) = sortedNumP(1);

index2 = 2;
index3 = 1;
while (index2 <= size(wSList,2)) && (index3 < numTrajs)
    % grab the next best guess
    tempW = sortedWS(:,index2);
    
    % check that it's not within the check radius of all better guesses
    ind = 1:index3;
    check1 = (tempW(1) < wSGuesses(:,ind)+searchBnds(1)) & (tempW(1) > wSGuesses(:,ind)-searchBnds(1));
    check2 = (tempW(2) < wSGuesses(:,ind)+searchBnds(2)) & (tempW(2) > wSGuesses(:,ind)-searchBnds(2));
    check3 = (tempW(3) < wSGuesses(:,ind)+searchBnds(3)) & (tempW(3) > wSGuesses(:,ind)-searchBnds(3));
    check = check1 & check2 & check3;
    
%     if length(sum(check)) > 1 || isempty(sortedNumP)
%         keyboard
%     end
    if (sum(sum(check)) > 0) && (sortedNumP(index2) > 1)
        wSGuesses(:,index3+1) = tempW;
        accNum(index3+1) = sortedNumP(index2);
        index3 = index3+1;
    end
    index2 = index2+1;
end

% remove unused guesses
index4 = accNum==0;
accNum(index4) = [];
wSGuesses(:,index4) = [];

%%% DEBUG %%%
% figure()
% hold on
% for ii = 1:size(wSList,2)
%     plot3(wSList(1,ii),wSList(2,ii),wSList(3,ii),...
%         'Marker','o','MarkerFaceColor','k')
% end
% view([1,1,1])
% hold off
% keyboard
%%%%%%%%%%%%%
end