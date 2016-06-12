function [peakList0, peakList] = FindMatchingPeaks(fitData0, fitData)
% FindMatchingPeaks: 

n0 = size(fitData0.hkl,1);

peakList0 = zeros(n0,1);
peakList = zeros(n0,1);
np = 1;

% for each peak in hkls0
for ii = 1:n0
    
    % find the matching hkl in hkls if it exists
    index1 = find(fitData.hkl(:,1) == fitData0.hkl(ii,1));
    index2 = find(fitData.hkl(:,2) == fitData0.hkl(ii,2));
    index3 = find(fitData.hkl(:,3) == fitData0.hkl(ii,3));
    
    index4 = intersect(intersect(index1,index2),index3);
    
    % if there is only one hit, record it
    if (length(index4)==1)
        peakList0(np) = ii;
        peakList(np) = index4;
        np = np+1;
        
    % if there are two hits, find the one with the closest omega value
    elseif length(index4)==2
        ome0 = fitData0.ome(ii);
        omes = fitData.ome(index4);
        deltaOme = abs(omes - ome0);
        [~,index5] = min(deltaOme);
        peakList0(np) = ii;
        peakList(np) = index4(index5);
        np = np+1;
    end
end

% remove trailing zeros
peakList0(peakList0==0)=[];
peakList(peakList==0) = [];

end