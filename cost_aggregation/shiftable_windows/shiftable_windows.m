% This approach aims at reducing the border localization
% problem of FW not constraining the support to be centered
% on the central position

% Example
% [D, C_min, C] = shiftable_windows(I1,I2,0,15,2,2,'SSD',1);

% Prepared by: Gabriel da Silva Vieira (Jan 2017)

function [D, C_min, C] = shiftable_windows(I1, I2, min_d, max_d, h, w, method, reverse)

% Execute block_matching to construct the DSI matrix
[D, C_min, C] = block_matching(I1 ,I2 ,min_d ,max_d, h, w, method, reverse);

C_forward = C;
C_backward = C;

for i=1:w
    shift_forward = imtranslate(C, [i 0]);
    shift_backward = imtranslate(C, [-i 0]);
    
    C_forward = min(C_forward, shift_forward);
    C_backward = min(C_backward, shift_backward);
end

    min_horizontal = min(C_forward, C_backward);
    C_down = min_horizontal;
    C_up = min_horizontal;

for i=1:h
    shift_down = imtranslate(min_horizontal, [0 i]); 
    shift_up = imtranslate(min_horizontal, [0 -i]);

    C_down = min(C_down, shift_down);
    C_up = min(C_up, shift_up);
end

    min_rectangle = min(C_down, C_up);

    [C_min, D] = min(min_rectangle, [], 3);
end
