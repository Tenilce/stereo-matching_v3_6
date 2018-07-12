% Block matching RGB, with raw cost, from image I1 to image I2,
% using integral images for computing the matching costs.

 % INPUT
 %   I1 the left stereo image
 %   I2 the right stereo image
 %   min_d minimum disparity
 %   max_d maximum disparity
 %   h, w heigth and width from the Fixed Windows, respectively 
 %   method used for calculating the correlation scores
 %   Valid values include: 'SAD', 'SSD', 'STAD', 'ZSAD', 'ZSSD', 'SSDNorm', 'NCC',
 %   'AFF', 'LIN', 'BTSSD', 'BTSAD'
 %   reverse used to calc disparity map from left to rigth, input 1 or -1,
 %   1 means regular disparity calculation, -1 means reverse disparity
 %   calculation
 %
 % OUTPUT
 %   D disparity values
 %   C_min cost associated with the minimum disparity at pixel (i,j)
 %   C  the cost volume for differences between I1 and I2
 %
  
 % Prepared by: Gabriel da Silva Vieira (Fev 2017)

function [D, C_min, C] = block_matching_RGB_raw_cost...
    (I1, I2, min_d, max_d, h, w, method, reverse)

[heigth, width, channels] = size(I1);
d_vals = min_d : max_d;
offsets = length(d_vals);
C_out = NaN(heigth, width, offsets, channels);

for i=1:channels
    [D1, C_min1, C1] = block_matching(I1(:,:,i), I2(:,:,i), min_d, max_d, h, w, method, reverse);
    C_out(:,:,:,i) = C1;
end

    C = min(C_out, [] ,4);
    [C_min, D] = min(C, [], 3);

end