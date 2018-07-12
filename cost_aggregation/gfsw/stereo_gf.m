% Guided Filter Support Weights - Stereo Vision

 % INPUT
 %   I1 the left stereo image
 %   I2 the right stereo image
 %   min_d minimum disparity
 %   max_d maximum disparity
 %   h, w heigth and width from the Fixed Windows, respectively 
 %   method used for calculating the correlation scores
 %   Valid values include: 'SAD', 'SSD', 'STAD', 'ZSAD', 'ZSSD', 'SSDNorm', 'NCC',
 %   'AFF', 'LIN', 'BTSSD', 'BTSAD', 'TAD_C+G'
 %   r - filter kernel has size r \times r
 %   eps - epsilon
 %   reverse used to calc disparity map from left to rigth, input 1 or -1,
 %   1 means regular disparity calculation, -1 means reverse disparity
 %   calculation
 %
 % OUTPUT
 %   D disparity values
 %   C_min cost associated with the minimum disparity at pixel (i,j)
 %   C  the cost volume for differences between I1 and I2
 %
 
 % Example:
 % [D, C_min, C] = stereo_gf(I1,I2,1,15,0,0,'TAD_C+G',19,0.0001,1);
 
 % Reference: 
 % C. Rhemann, A. Hosni, M. Bleyer, C. Rother, M. Gelautz, 
 % Fast Cost-Volume Filtering for Visual Correspondence and Beyond, CVPR11
 
 % Prepared by: Gabriel da Silva Vieira (Nov 2017)


function [D, C_min, C] = stereo_gf(I1, I2, min_d, max_d, h, w, method, r, eps, reverse)

% Execute block_matching to construct the DSI matrix
[D1, C_min1, C] = block_matching(I1 ,I2 ,min_d ,max_d, h, w, method, reverse);

% the range of disparity values from min_d to max_d inclusive
d_vals = min_d : max_d;
offsets = length(d_vals);

for off=1:offsets
    disp = C(:,:,off);
    disp_filtered = imguidedfilter(disp, double(I1)/255, ...
        'NeighborhoodSize', [r r], 'DegreeOfSmoothing', eps);
    
    C(:,:,off) = disp_filtered;
end

[C_min, D] = min(C, [], 3);

end

