% Function which works as IMAGESC(d)
%
% Input: disparity map
% Output: representation of the disparity map
%

function [d_color] = display_dmap(d)

% Map disparity into the range [0,1]
d_scaled = (d-min(d(:)))/range(d(:));
% Colorize occluded pixels to be blue
d_color = repmat(d_scaled(:),[1 3]);
occ_inds = isnan(d(:));
d_color(occ_inds,1) = 0;
d_color(occ_inds,2) = 0;
d_color(occ_inds,3) = 1;
d_color = reshape(d_color,[size(d) 3]);
% Display image
image(d_color);

end