function binaryImg = BuildBinaryImg(imgInfo, instr, options)
% Used for debugging a section of TotalSearchGrainFinderLR
% 


% defaults
threshold = 1000;

if exist('options','var')
    if isfield(options,'threshold')
        threshold = options.threshold;
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
end