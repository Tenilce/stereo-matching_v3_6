% Perform Block Matching with Dynamic Programming.
%  Adapted from the example on the Matlab help site.
%clear all
close all

% Read images.
leftI = rgb2gray(imread('scene1.row3.col1.ppm'));
rightI = rgb2gray(imread('scene1.row3.col2.ppm'));
figure(1), imshow(leftI, []), title('Left image');
figure(2), imshow(rightI, []), title('Right image');

leftI = single(leftI);
rightI = single(rightI);
nRowsLeft = size(leftI, 1);

% Blocksize is (2*halfBlockSize+1)-by-(2*halfBlockSize+1).
halfBlockSize = 3;
blockSize = 2*halfBlockSize+1;

disparityRange = 15;


% Basic block matching chooses the optimal disparity for each pixel based
% on its own cost function alone. Now we want to allow a pixel to have a
% disparity with possibly sub-optimal cost for it locally. This extra cost
% must be offset by increasing that pixel's agreement in disparity with its
% neighbors. In particular, we constrain each disparity estimate to lie
% with $\pm 3$ values of its neighbors' disparities, where its neighbors
% are the adjacent pixels along an image row. The problem of finding the
% optimal disparity estimates for a row of pixels now becomes one of
% finding the "optimal path" from one side of the image to the other. To
% find this optimal path, we use the underlying block matching metric as
% the cost function and constrain the disparities to only change by a
% certain amount between adjacent pixels. This is a problem that can be
% solved efficiently using the technique of dynamic programming [3,4].

D = zeros(size(leftI), 'single');
finf = 1e3; % False infinity
disparityCost = finf*ones(size(leftI,2), 2*disparityRange + 1, 'single');
disparityPenalty = 0.5; % Penalty for disparity disagreement between pixels
hWaitBar = waitbar(0,'Using dynamic programming for smoothing...');
% Scan over all rows.
for m=1:nRowsLeft
    disparityCost(:) = finf;
    % Set min/max row bounds for image block.
    minr = max(1,m-halfBlockSize);
    maxr = min(nRowsLeft,m+halfBlockSize);
    % Scan over all columns.
    for n=1:size(leftI,2)
        minc = max(1,n-halfBlockSize);
        maxc = min(size(leftI,2),n+halfBlockSize);
        % Compute disparity bounds.
        mind = max( -disparityRange, 1-minc );
        maxd = min( disparityRange, size(leftI,2)-maxc );
        % Compute and save all matching costs.
        for d=mind:maxd
            disparityCost(n, d + disparityRange + 1) = ...
                sum(sum(abs(leftI(minr:maxr,(minc:maxc)+d) ...
                - rightI(minr:maxr,minc:maxc))));
        end
    end
    
    % Process scan line disparity costs with dynamic programming.
    optimalIndices = zeros(size(disparityCost), 'single');
    cp = disparityCost(end,:);
    for j=size(disparityCost,1)-1:-1:1
        % False infinity for this level
        cfinf = (size(disparityCost,1) - j + 1)*finf;
        % Construct matrix for finding optimal move for each column
        % individually.
        [v,ix] = min([cfinf cfinf cp(1:end-4)+3*disparityPenalty;
                      cfinf cp(1:end-3)+2*disparityPenalty;
                      cp(1:end-2)+disparityPenalty;
                      cp(2:end-1);
                      cp(3:end)+disparityPenalty;
                      cp(4:end)+2*disparityPenalty cfinf;
                      cp(5:end)+3*disparityPenalty cfinf cfinf],[],1);
        cp = [cfinf disparityCost(j,2:end-1)+v cfinf];
        % Record optimal routes.
        optimalIndices(j,2:end-1) = (2:size(disparityCost,2)-1) + (ix - 4);
    end
    % Recover optimal route.
    [~,ix] = min(cp);
    D(m,1) = ix;
    for k=1:size(D,2)-1
        D(m,k+1) = optimalIndices(k, ...
            max(1, min(size(optimalIndices,2), round(D(m,k)) ) ) );
    end
    waitbar(m/nRowsLeft, hWaitBar);
end
close(hWaitBar);
D = D - disparityRange - 1;


% Show resulting disparity map in pseudocolor.
% For display purposes, we saturate the depth map to have only positive
% values. In general, slight angular misalignment of the stereo cameras
% used for image acquisition can allow both positive and negative
% disparities to appear validly in the depth map. In this case, however,
% the stereo cameras were near perfectly parallel, so the true disparities
% have only one sign. Thus this correction is valid.
figure
imshow(D,[]), axis image, colormap('gray'), colorbar, impixelinfo;
caxis([0 disparityRange]);
title('Disparity map from block matching with dynamic programming');
