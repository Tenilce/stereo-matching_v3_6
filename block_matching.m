% Block matching from image I1 to image I2,
% using integral images for computing the matching costs.

 % INPUT
 %   I1 the left stereo image
 %   I2 the right stereo image
 %   min_d minimum disparity
 %   max_d maximum disparity
 %   h, w heigth and width from the Fixed Windows, respectively 
 %   method used for calculating the correlation scores
 %   Valid values include: 'SAD', 'SSD', 'STAD', 'ZSAD', 'ZSSD', 'SSDNorm', 'NCC',
 %   'AFF', 'LIN', 'BTSSD', 'BTSAD', 'TAD_C+G'
 %   reverse used to calc disparity map from left to rigth, input 1 or -1,
 %   1 means regular disparity calculation, -1 means reverse disparity
 %   calculation
 %
 % OUTPUT
 %   D disparity values
 %   C_min cost associated with the minimum disparity at pixel (i,j)
 %   C  the cost volume for differences between I1 and I2
 %
 
 % References: 
 % Gabriele Facciolo, Nicolas Limare, Enric Meinhardt. Integral Images for Block Matching, 2014
 
 % Example
 % [D, C_min, C] = block_matching(I1,I2,0,15,2,2,'SAD',1);
 % using STAD
 % [D, C_min, C] = block_matching(I1,I2,0,15,2,2,['STAD', {20}],1);
 
 % Prepared by: Gabriel da Silva Vieira (Jan 2017)

function [D, C_min, C] = block_matching(I1, I2, min_d, max_d, h, w, method, reverse) 

% Prepared to use 'STAD'
if size(method,2) == 2
    truncated_value = method{2};
    method = method(1);
end

I1 = double(I1) / 255;
I2 = double(I2) / 255;

% window size to aggerate the cost
w_heigth = h*2+1;
w_width = w*2+1;

[h, w, channels] = size(I1);

% the range of disparity values from min_d to max_d inclusive
d_vals = min_d : max_d;
offsets = length(d_vals);
C = ones(h,w,offsets); % the cost volume

% validate input arguments
[I1, I2] = valid_inputs(I1, I2, offsets, w_heigth, w_width, method, reverse);

% Precomputed images needed for the cost evaluation (μ, σ,...)
if strcmp(method,'ZSSD') || strcmp(method,'SSDNorm') || strcmp(method,'NCC')...
        || strcmp(method,'AFF') || strcmp(method,'LIN')
    meanI1 = NaN(h,w,channels);
    normI1 = NaN(h,w,channels);
    varI1 = NaN(h,w,channels);

    meanI2 = NaN(h,w,channels);
    normI2 = NaN(h,w,channels);
    varI2 = NaN(h,w,channels);

    for k=1:channels
        [mI1, nI1, vI1] = stuff(I1(:,:,k), w_heigth, w_width);
        meanI1(:,:,k) = mI1;
        normI1(:,:,k) = nI1;
        varI1(:,:,k) = vI1;

        [mI2, nI2, vI2] = stuff(I2(:,:,k), w_heigth, w_width);
        meanI2(:,:,k) = mI2;
        normI2(:,:,k) = nI2;
        varI2(:,:,k) = vI2;
    end
end

% Prepare same constante values
if strcmp(method, 'STAD') 
    if isempty(truncated_value)
        error('A truncated value needs to be informed');
    end
    %truncated_value = input('!!!!!! >> Inform a value to be truncated: ');
    % p.e, threshold = 20, h=w=7, and compare with SAD with h=w=7 
elseif strcmp(method, 'ZSAD') % Zero-mean SAD
     I1 = I1 - meanI1;
     I2 = I2 - meanI2;
elseif strcmp(method, 'BTSSD') || strcmp(method, 'BTSAD')
    I1_neg = (I1 + (imtranslate(I1, [-1 0])))./2;
    I1_pos = (I1 + (imtranslate(I1, [1 0])))./2;
    [I1_max, I1_min] = max_min_3(I1_neg, I1_pos, I1);
elseif strcmp(method, 'TAD_C+G')
    thresColor = 7/255;     
    thresGrad = 2/255;      
    gamma = 0.11;           % (1- \alpha) 
    
    if channels > 1
        I1_ = sum(I1,channels) / channels; 
        I2_ = sum(I2,channels) / channels; 
    else
        I1_ = I1;
        I2_ = I2;
    end
        
    fx_l = gradient(I1_);
    fx_l = fx_l + 0.5; % To get a range of values between 0 to 1
        
    fx_r = gradient(I2_);
    fx_r = fx_r + 0.5; % To get a range of values between 0 to 1

end

% To use sum function in an appropriate way
if channels == 1
    dimension = 3;
else
    dimension = channels;
end

% the main loop
for off=1:offsets    
    d = d_vals(off);
    v_shift = imtranslate(I2, [(reverse)*d 0]);
    
    % precompute pixel distances
    if strcmp(method, 'SAD') || strcmp(method, 'ZSAD')
        dist = abs(I1 - v_shift); % Sum of Absolute Differences (SAD)
        dist = sum(dist,dimension)/channels;
        C(:,:,off) = sum_patches(dist,w_heigth,w_width);
    elseif strcmp(method, 'SSD')
        dist = (I1 - v_shift).^2; % Sum of Squared Differences (SSD)
        dist = sum(dist,dimension)/channels;
        C(:,:,off) = sum_patches(dist,w_heigth,w_width);
    elseif strcmp(method, 'STAD')
        dist = abs(I1 - v_shift); % Sum of Absolute Differences (SAD)
        dist = sum(dist,dimension)/channels;
        dist(dist > truncated_value) = truncated_value; % Sum of truncated absolute differences (STAD)
        C(:,:,off) = sum_patches(dist,w_heigth,w_width);
    elseif strcmp(method, 'ZSSD')
        dist = (I1 - v_shift).^2;
        dist = sum(dist,dimension)/channels;
        ssd = sum_patches(dist,w_heigth,w_width);
        C(:,:,off) = ssd - (meanI1 - meanI2).^2;
    elseif strcmp(method, 'SSDNorm')
        dist = I1.*v_shift;
        
        prod = NaN(h,w,channels);
        norm_v_shift = NaN(h,w,channels);
        
        for k=1:channels
            prod(:,:,k) = sum_patches(dist(:,:,k), w_heigth, w_width);
            norm_v_shift(:,:,k) = norm_rectangle(v_shift(:,:,k), w_heigth, w_width);
        end

        den = normI1.*norm_v_shift;
        % prevent divisions by zero, Under that circumstance the cost should be set to 2
        den(den==0) = 1; prod(den==0) = 0; 
        result = 2 - (2.*(prod./den));
        C(:,:,off) = sum(result,dimension)/channels; % SSD/Norm
    elseif strcmp(method, 'NCC')
        dist = I1.*v_shift;
        
        prod = NaN(h,w,channels);
        mean_v_shift = NaN(h,w,channels);
        var_v_shift = NaN(h,w,channels);
        
        for k=1:channels
            prod(:,:,k) = (sum_patches(dist(:,:,k), w_heigth, w_width)) / (w_heigth * w_width);
            mean_v_shift(:,:,k) = mean_rectangle(v_shift(:,:,k), w_heigth, w_width);
            var_v_shift(:,:,k) = variance_rectangle(v_shift(:,:,k), w_heigth, w_width);
        end
        
        medias = meanI1 .* mean_v_shift;
        num = prod - medias;
        
        den = sqrt(varI1 .* var_v_shift);
        % prevent divisions by zero, Under that circumstance the cost should be set to 1
        %num(num==0) = 1; den(den==0) = 1;
        
        corr = num ./ den; % Normalized Cross Correlation (NCC)
        corr(isnan(corr)) = -1;
        
        result = sum(corr,dimension)/channels;
        
        C(:,:,off) = 1 - result;
    elseif strcmp(method, 'AFF')
        dist = I1.*v_shift;
        
        for k=1:channels
            prod = (sum_patches(dist(:,:,k), w_heigth, w_width)) / (w_heigth * w_width);
            mean_v_shift = mean_rectangle(v_shift(:,:,k), w_heigth, w_width);
            var_v_shift = variance_rectangle(v_shift(:,:,k), w_heigth, w_width);
            
            medias = meanI1(:,:,k) .* mean_v_shift;
            num = prod - medias;
            den = sqrt(varI1(:,:,k) .* var_v_shift);
            corr = num ./ den;
            
            result = max(varI1(:,:,k),var_v_shift).*min(1,1-(corr.*abs(corr))); % "affine" similarity measure          
        end
        
        C(:,:,off) = result;
    elseif strcmp(method,'LIN')
        dist = I1.*v_shift;
        
        result = NaN(h,w,channels);
        
        for k=1:channels
            prod = (sum_patches(dist(:,:,k), w_heigth, w_width)) .^ 2;
            norm_v_shift = norm_rectangle(v_shift(:,:,k), w_heigth, w_width);
            
            den = (normI1(:,:,k).*norm_v_shift).^2;
            den(den==0) = 1; 
            
            f1 = (max(normI1(:,:,k),norm_v_shift)).^2;
            f2 = 1 - (prod./den);
            
            result(:,:,k) = f1.*f2;  % This is a simpler variant of the AFF cost         
        end

        C(:,:,off) = sum(result,dimension) / channels;
    elseif strcmp(method, 'BTSSD')
        I2_neg = (v_shift + (imtranslate(v_shift, [-1 0])))./2;
        I2_pos = (v_shift + (imtranslate(v_shift, [1 0])))./2;
        [I2_max, I2_min] = max_min_3(I2_neg, I2_pos, v_shift);
        
        result = NaN(h,w,channels);
        for k=1:channels      
            d_LR = max_min_3(zeros(h,w), (I1(:,:,k) - I2_max(:,:,k)), (I2_min(:,:,k) - I1(:,:,k)));
            d_RL = max_min_3(zeros(h,w), (v_shift(:,:,k) - I1_max(:,:,k)), (I1_min(:,:,k) - v_shift(:,:,k)));      
            d_bt = (min(d_LR,d_RL)).^2; % Birchfield and Tommasi Sampling Insensitive distance BTSSD
            result(:,:,k) = sum_patches(d_bt, w_heigth,w_width);
        end
            
        C(:,:,off) = sum(result,dimension) / channels;
    elseif strcmp(method, 'BTSAD')
        I2_neg = (v_shift + (imtranslate(v_shift, [-1 0])))./2;
        I2_pos = (v_shift + (imtranslate(v_shift, [1 0])))./2;
        [I2_max, I2_min] = max_min_3(I2_neg, I2_pos, v_shift);
        
        result = NaN(h,w,channels);
        for k=1:channels      
            d_LR = max_min_3(zeros(h,w), (I1(:,:,k) - I2_max(:,:,k)), (I2_min(:,:,k) - I1(:,:,k)));
            d_RL = max_min_3(zeros(h,w), (v_shift(:,:,k) - I1_max(:,:,k)), (I1_min(:,:,k) - v_shift(:,:,k)));      
            d_bt = abs(min(d_LR,d_RL)); % Birchfield and Tommasi Sampling Insensitive distance BTSAD
            result(:,:,k) = sum_patches(d_bt, w_heigth,w_width);
        end
            
        C(:,:,off) = sum(result,dimension) / channels;
    elseif strcmp(method, 'TAD_C+G')
       
        % Truncated SAD of color images for current displacement
        p_color = abs(I1 - v_shift); % Absolute Differences (AD)
        p_color = sum(p_color,dimension) / channels;
        p_color = min(p_color,thresColor);
        
        tmp = imtranslate(fx_r,[(reverse)*d 0]);
        p_grad = abs(tmp - fx_l);
        p_grad = min(p_grad,thresGrad);
    
        p = gamma*p_color + (1-gamma)*p_grad; % Combined color and gradient
        
        C(:,:,off) = sum_patches(p,w_heigth,w_width);
    end

end
    
    [C_min, D] = min(C, [], 3);        

end

% Function to determine max and min values between 3 matrices
%
function [I_max, I_min] = max_min_3(I1,I2,I3)

 [h, w, c] = size(I1);
 I_max = NaN(h,w,c);
 I_min = NaN(h,w,c);
     for k=1:c
         I_max(:,:,c) = max(max(I1(:,:,c),I2(:,:,c)),I3(:,:,c));
         I_min(:,:,c) = min(min(I1(:,:,c),I2(:,:,c)),I3(:,:,c));
     end
end


function [I1, I2] = valid_inputs(I1, I2, offsets, r1, r2, method, reverse)

[h1,w1] = size(I1);
[h2,w2] = size(I2);

% Check to see if both the left and right images have same number of rows
% and columns
if h1~=h2 || w1~=w2
    error('Both left and right images should have the same number of rows and columns');
end

if offsets < 1
    error('Offsets need to be iqual or greater than 1');
end

if r1<1 || r2<1
    error('r1 and r2 need to be iqual or greater than 1');
end

if isequal(reverse,-1)
    I = I1;
    I1 = I2;
    I2 = I;
end

end
