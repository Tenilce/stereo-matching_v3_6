function D = cc_min_disps(I1, D1, D2)


D1(isnan(D1)) = 1000;
D2(isnan(D2)) = 500;

[heigth, width, ~] = size(I1);
Y = repmat((1:heigth)', [1 width]);
X = repmat(1:width, [heigth 1]) - D1;
X(X<1) = 1;
indices = sub2ind([heigth,width],Y,X);

D = nanmin(D1,D2(indices));

end