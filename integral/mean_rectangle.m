% Function to compute mean value of a patch
% Given an image, the sum of value of pixels in the rectangular region
% dividing by the quantity of number of pixels gives us the mean value of pixels in the region.

% Prepared by: Gabriel da Silva Vieira (Jan 2017)

function m = mean_rectangle(I, h, w)

I = make_border(I,h,w);

ii = integral_image(I);

[h_ii, w_ii] = size(ii);
rec = NaN(h_ii, w_ii);

%prepare bordes 
p_h = (h-1)/2;
q_h = p_h+1;
p_w = (w-1)/2;
q_w = p_w+1;

%Start loops
begin_h = q_h+1;
begin_w = q_w+1;
max_h = h_ii - h;
max_w = w_ii - w;

r1_floor = floor(h/2);
r2_floor = floor(w/2);

for i=begin_h:max_h
    for j=begin_w:max_w
        rec(i+r1_floor,j+r2_floor) = sum_rectangle(ii, i, j, h, w);
    end
end

m = rec/(h*w);

m = crop_border(m,h,w);

end