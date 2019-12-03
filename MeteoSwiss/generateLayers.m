function W_multi = generateLayers(stations)
% GENERATELAYERS genrates 2-binary layers wrt gps and altitude proximity
% for the given set of nodes on stations.st_indx
% W_multi : cell output of the Adjacency matrices of the layers

sparse = 0.1;
N = length(stations.st_indx);
num_edge = round(sparse*N*(N-1)); % x2
gps = stations.GPS(stations.st_indx,:);
alt = stations.altitude(stations.st_indx,:);

%% GPS
Dist = pdist2(gps,gps,'euclidean');
Dist(logical(eye(size(Dist)))) = max(Dist(:));
gps_dist = sort(Dist(:));
thresh_gps = gps_dist(num_edge);
W1 = double(Dist <= thresh_gps);

%% Altitude
Dist = pdist2(alt,alt,'euclidean');
Dist(logical(eye(size(Dist)))) = max(Dist(:));
alt_dist = sort(Dist(:));
thresh_alt = alt_dist(num_edge);
W2 = double(Dist <= thresh_alt);

%% Normalize the volume wrt edgeset
W1 = W1*N/num_edge;
W2 = W2*N/num_edge;

W_multi = {W1,W2};