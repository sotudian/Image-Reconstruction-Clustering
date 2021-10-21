function w = FindWeights (X, U, V, m, K)
%Finds weights for PCM Clustering.
c = size (V, 1);
w = zeros (c, 1);
mf = U.^m;
dist = Distance_Function (V, X);
dist = dist .^ 2;
w = sum((mf.*dist)') ./ sum(mf');
