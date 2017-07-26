function cleanData = cleanData(data,distance)
% CLEANDATA This function cleans the data from a point cloud matrix, which
% should be (number of points) x 6, where the first three columns are
% coordinates and the last three columns are colors. 'Cleaning' the data
% means removing data containing Inf and NaN values and applying a distance
% threshold.

% Remove Inf's and NaN's
data = data(~any(isnan(data),2),:);
data = data(~any(isinf(data),2),:);
% Remove points that are far away
dists = sqrt(data(:,1).^2 + data(:,2).^2 + data(:,3).^2);
cleanData = data(dists <= distance,:);