function grains = TotalSearchGrainFinderLR2(imgInfo, n, omeBnds, omeStep, instr, material, options)
% TotalSearchGrainFinderLR2: Finds grain orientations by the same method as
% 'TotalSearchGrainFinder', but sacrifices speed to reduce the required
% RAM.  Adds parallel processing compatability.
% 
% USAGE: 
%   grains = TotalSearchGrainFinderLR2(imgInfo, n, omeBnds, omeStep,...
%       instr, material, options);
%
% AUHTORS: Timothy Long, Darren Pagan*, and Mark Obstalecki*
%
% INPUTS:
%   imgInfo is a structure:
%   Contains the information needed to load the images.
%       .dir is a string 
%           complete path to the fast sweep file, e.x. 'C:/data/'
%
%       .stem is a string
%           stem of the file, e.x. 'FF_' for far field images
%
%       .sweepNumArray is m x 1
%           vector containing the image numbers of the fast sweep files
%
%       .numDark is 1 x 1
%           number of dark images at the front of the fast sweep file       
% 
%       .numPerStack is 1 x 1
%           numbers of data images in the fast sweep file
%
%   n is a natural number:
%       Controls the density of the grid over orientation space.  A larger
%       value of n increases the number of points in the grid.  Must be
%       a positive number!  Non-interger values may give highly uneven grid
%       spacing.
%
%   omeBnds is 1 x 2
%         the omega bounds of the diffraction experiment
%
%   omeStep is 1 x 1
%         the omega integrating step size of the diffractometer
%   
%   instr is a structure:
%       Contains fields describing the insturmentation used for the
%       experiment.
%
%   material is a structure:
%       Contains fields describing the material being studied in the
%       experiment.  Non-standard fields are:
%
%       .sym is s (# of symmetry operators) x 4:
%           The symmetry operators of the crystal lattice, represented as
%           quaterions.  If this field is not provided, only the identity
%           operator (1,0,0,0) will be used.
%
%   options is a structure:
%       Contains fields that change parameters or alter code execution.
%       Possible fields are:
%
%       .threshold is a scalar:
%           Used to convert the detector images to binary images.  Pixels
%           with measured intensities above this value are 1, all other
%           pixels are 0.  Default is 1000.
%
%       .searchAdjacentPix is a boolean:
%           If true, pixels adjacent to the pixel the peak is predicted to
%           fall on will be checked for intensity during the primary search
%           step.  Default is false.
%
%       .completenessTol is a scalar:
%           The minimum completeness a test orientation must have to become
%           a candidate orientation.  Must be in the range (0,1).  Default
%           is 0.75 (75%).
%
%       .searchRange is a scalar:
%           The multiple of the maximum grid spacing (2/n) away from each
%           candidate orientation that the post-processing looks for higher
%           completion candidate orientations when searching for grains.
%           Values less than 1 may cause unusual results.  Default is 1.1
%           (should be nearest neighbors only).
% 
%       .numSimRings is a scalar:
%           The number of rings to look for peaks on.
%
% NOTES:
%   *T.L. copied code by D.P. and M.O. in writing this function.  T.L. is
%   solely responsible for any bugs, errors, etc in this function.
%
%   Currently only supports .ge2 image files


% Internal Notes:
%   The binary images are stored as a cell array of sparse matricies
%   because matlab only has 1 or 2-D sparse arrays.


% check input
if n <= 0
    error('Grid refinment parameter (n) must be a positive number.')
end


% defaults
threshold = 1000;
searchAdjacentPix = false;
completenessTol = 0.75;
searchRange = 1.1;
hkl = GethklsTruncate(instr,material,material.tth(end));
hkl = RemoveZeroIntensityHKL(material,hkl);

sym = [1,0,0,0];

% overwrite defaults if passed in
if isfield(material,'sym')
    sym = material.sym;
end

if exist('options','var')
    if isfield(options,'threshold')
        threshold = options.threshold;
    end
    
    if isfield(options,'searchAdjacentPix')
        searchAdjacentPix = options.searchAdjacentPix;
    end
    
    if isfield(options,'completenessTol')
        completenessTol = options.completenessTol;
    end
    
    if isfield(options,'searchRange')
        searchRange = options.searchRange;
    end
    
    if(isfield(options,'numSimRings'))
       numSimRings=options.numSimRings;
       hkl = GethklsTruncate(instr,material,material.tth(numSimRings));
       hkl = RemoveZeroIntensityHKL(material,hkl);
    end
end


%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
% build binary images
disp('Converting images to binary...')

numImg = size(imgInfo.sweepNumArray,1);
nPS = imgInfo.numPerStack;

binaryImg = mat2cell(zeros(numImg*nPS,1),numImg*nPS,1);
for ii = 1:numImg
    disp(['Processing fastsweep #' num2str(ii) ' of ' num2str(numImg)])
    % load one fastsweep
    imageStack = FastsweepRamLoad(imgInfo.dir, imgInfo.stem,...
        imgInfo.sweepNumArray(ii), imgInfo.numDark,...
        nPS, instr, options, []);
    
    % Convert each image to a sparse binaries and store it
    for jj = 1:size(imageStack,3)
        binaryImg{nPS*(ii-1) + jj} = ...
            sparse(imageStack(:,:,jj)>threshold);
    end
end


%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
% Search all points in the fundamental region
disp('Starting primary search...')
gridSize = (n+1)^3*8;

% dist = 2/n;
% list = (-1:dist:1);

% rowNum = 0;
% colNum = 1;
% sheetNum = 1;
% faceNum = 1;

index1 = 1;
hits = struct();
for ii=1:gridSize
    inFndReg = false;
    
    % new grid generator
    q = GenerateQuatGridLR(n,ii);

    % project the point onto the unit 4-sphere
    q = q/norm(q);

    % check if the point is in the fundamental region
    if size(sym,1)>1
        delta = [sym(:,1)-q(1), sym(:,2)-q(2), sym(:,3)-q(3), sym(:,4)-q(4)];
        dist = sum(delta.^2,2).^0.5;
        [~,zoneNum] = min(dist);

        if zoneNum==1
            inFndReg = true;
        end
    else
        inFndReg = true;
    end

    % only search points that are in the fundamental region
    if inFndReg
        Rsc = RMatOfQuat(q');

        % run virtual diffractometer on grain
        dataTemp=RotDiffractionExp(instr.ki,hkl,Rsc,material.b,omeBnds);
        % find frame numbers for each peak
        [frameno,~] = FindFrameNo(dataTemp.ome,omeBnds,omeStep);
        % pixel row and column on the detector
        [~,pixCoords] =...
            DetDetectorIntercept(instr,dataTemp.koL',zeros(size(dataTemp.koL')));

        % remove peaks that aren't on the detector
        index=(pixCoords(:,2)>(1) & ...
            pixCoords(:,2)<(instr.det.pixelDim(1))  & ...
            pixCoords(:,1)>(1)  & ...
            pixCoords(:,1)<(instr.det.pixelDim(2)));
        pixCoords = pixCoords(index,:);
        frameno = frameno(index);

        % check for intensity at those pixels
        pixVals = zeros(size(frameno));
        
        index2 = sub2ind([instr.det.pixelDim(2), instr.det.pixelDim(1)],...
            pixCoords(:,2), pixCoords(:,1));
        
        % check adjacent pixels
        if searchAdjacentPix
            i2a = sub2ind([instr.det.pixelDim(1), instr.det.pixelDim(2)],...
                min(pixCoords(:,1)+1,instr.det.pixelDim(1)), pixCoords(:,2));

            i2b = sub2ind([instr.det.pixelDim(1), instr.det.pixelDim(2)],...
                max(pixCoords(:,1)-1,1), pixCoords(:,2));

            i2c = sub2ind([instr.det.pixelDim(1), instr.det.pixelDim(2)],...
                pixCoords(:,1), min(pixCoords(:,2)+1,instr.det.pixelDim(2)));

            i2d = sub2ind([instr.det.pixelDim(1), instr.det.pixelDim(2)],...
                pixCoords(:,1), min(pixCoords(:,2)-1,1));

            for jj = 1:length(frameno)
                pixVals = [binaryImg{frameno(jj)}(index2(jj)),...
                    binaryImg{frameno(jj)}(i2a(jj)),...
                    binaryImg{frameno(jj)}(i2b(jj)),...
                    binaryImg{frameno(jj)}(i2c(jj)),...
                    binaryImg{frameno(jj)}(i2d(jj))];
                pixVals = max(pixVals,[],2);
            end
        else
            for jj = 1:length(frameno)
                pixVals(jj) = binaryImg{frameno(jj)}(pixCoords(jj,2),pixCoords(jj,1));
            end
        end

        % calculate completeness
        comp = sum(pixVals)/length(pixVals);
        
        %%% DEBUG %%%
%         disp(['Predicted ' num2str(length(frameno)) ' peaks, found '...
%             num2str(comp*length(frameno))])
%         disp(['Orientation: ' num2str(q)])
%         disp('-------')
%         disp(' ')
        %%%%%%%%%%%%%

        % record hits
        if comp > completenessTol
            %%% DEBUG %%%
            disp(['Hit #' num2str(index1)])
            disp(['Orientation = ' num2str(q)])
            disp(['Completeness = ' num2str(comp)])
            disp('-----')
            disp(' ')
            
            hits(index1).q = q';
            hits(index1).comp = comp;
            index1 = index1+1;
        end
    end
end

% remove duplicate hits
%   Duplicate hits are possible because I include all points in each cubic
%   cell of the tesseract, so all verticies, all edges, and some faces are
%   repeated in the array of test points.?
pts = [hits.q];
[~,index3] = unique(pts','rows');
hits = hits(index3);


disp('Primary processing is complete')
disp(['Found ' num2str(size(hits,2)) ' candidate orientations'])
disp(' ')


%+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
% find local maxima among the hits
disp('Starting post processing')
disp(' ')

% needed variables
dist0 = 2/n * searchRange;
pts = [hits.q];
index4 = 1;

grains = struct(); % create storage variable
for ii = 1:size(hits,2)
    testQ = hits(ii).q;
    deltas = [pts(1,:) - testQ(1);
              pts(2,:) - testQ(2);
              pts(3,:) - testQ(3);
              pts(4,:) - testQ(4)];
    
    % find nearyby hits
    dist = sum(deltas.^2,1).^0.5;
    index5 = dist < dist0;
    
    % Check if the point has the highest completeness among its
    % neighbors.  If it does, record it.
    if ~isempty(index5) && (sum([hits(index5).comp] > hits(ii).comp)) == 0
        grains(index4).Rsc = RMatOfQuat(testQ);
        grains(index4).r = RodOfQuat(testQ);
        grains(index4).completeness = hits(ii).comp;
        index4 = index4 + 1;
    end
end

% final status
disp('Completed indexing.')
disp(['Found ' num2str(size(grains,2)) ' grains.'])
disp('Parameters were: ')
disp(['Intensity threshold: ' num2str(threshold)])
disp(['Completeness tolerance:' num2str(completenessTol)])
disp(['Square grid spacing: ' num2str(2/n)])
disp('-----')
disp(' ')
end