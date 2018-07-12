% Left-Right Consistency Check

 % INPUT
 %   I1 the left stereo image - Reference Image
 %   min_d minimum disparity
 %   max_d maximum disparity
 %   D1 initial disparity map from left to right
 %   D2 initial disparity map from rigth to left
 %   gamma_c parameter to adjust spatial and color similarity in eq. (6)
 %   gamma_d parameter to adjust spatial and color similarity in eq. (6)
 %   r_median filter dimension of weighted median, has size r_median \times
 %   r_median in eq. (6)
 %
 % OUTPUT
 %   D_filled disparity values filled with a weighted bilateral filter
 %   D_occ disparity map which shows occused points
 %
 
 % Reference: 
 % C. Rhemann, A. Hosni, M. Bleyer, C. Rother, M. Gelautz, 
 % Fast Cost-Volume Filtering for Visual Correspondence and Beyond, CVPR11
 
 % Example
 % [D_filled, D_occ] = cc_min(I1,D1,D2);
 
 % Prepared by: Gabriel da Silva Vieira (Nov 2017)

function [D_filled, D_occ] = cc_min(I1, D1, D2)

D1(D1<=0) = 1;
D2(D2<=0) = 1;
D1(isnan(D1)) = 1000;
D2(isnan(D2)) = 500;

[h, w, ~] = size(I1);

% Left-right consistency check
Y = repmat((1:h)', [1 w]);
X = repmat(1:w, [h 1]) - D1;
X(X<1) = 1;
indices = sub2ind([h,w],Y,X);

D_occ = D1;
D_occ(abs(D1 - D2(indices))>=1) = -1;

% Fill and filter (post-process) pixels that fail the consistency check
D_filled = fill_min(D_occ);

end