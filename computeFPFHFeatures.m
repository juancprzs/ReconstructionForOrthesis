function fpfh = computeFPFHFeatures(coords,normals,radius,bins)
% COMPUTEFPFHFEATURES This function computes the Fast Point Feature
% Histograms as described in Fast Point Feature Histograms (FPFH) for 3D 
% Registration, by Bogdan et al. The parameters are defined as follows:
% - COORDS is a N-by-3 matrix of the coordinates of the points in the point 
% cloud as registered by a Kinect sensor.
% - NORMALS is a N-by-3 matrix representing the normal vectors of each of
% the points given in the COORDS matrix. Naturally, they must correspond in
% the order.
% - RADIUS scalar value representing the radius in which neighbors of a
% point will be searched. The time this function takes is extremely
% sensible to this parameter, so this must be kept in mind when thinking
% about computational resources.
% - BINS number of bins per histogram. Keep in mind that three features are
% computed and then a histogram is computed, therefore the actual number of
% bins of the resulting histogram will be 3*bins.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compute simplified features (SPFH)
% Preallocate memory for the features
fts = zeros(size(coords,1),3*bins);
% Perform range search to find all the neighbors (answer comes as cell
% array). It takes about one minute!
fprintf('Beginning to look for neighbors. May take a while!\n...')
begin = tic;
[neighbors,dist] = rangesearch(coords,coords,radius);
fprintf('Done looking for neighbors. Total of %g seconds\n',toc(begin))

% Problem here: all points will find themselves as close neighbors, and
% they will appear as the first entry of each array in the cell array
% 'neighbors', so we should avoid that by making iterations only from the
% second entry of each of the arrays within 'neighbors'.
% Note: I tried deleting the entry with 'cellfun(@(x)x(2:end),neighbors)'
% which should be fast, but it is NOT.

% Define 'edges' for the histograms, based in the number of bins 
alphaEdges  = linspace(-1,1,bins+1);
phiEdges    = alphaEdges;
thetaEdges  = linspace(-pi,pi,bins+1);

% Iterate through ALL the points (about 24 minutes) with for-loop
fprintf('Beginning to compute SPFH\n')
t1 = tic;
parfor idx = 1:size(coords,1)
    % Normal vector of the current point
    currentNormalVec = normals(idx,:);
    % Coordinates of the current point
    currentCoordinate = coords(idx,:);
    
    % Extract neighbors according to rangesearch
    currentNeighbors = neighbors{idx};
    
    % Save some memory for allocating the features (and so that computing
    % the histogram becomes easier!)
    alphas  = zeros(size(currentNeighbors,2)-1,1);
    phis    = zeros(size(currentNeighbors,2)-1,1);
    thetas  = zeros(size(currentNeighbors,2)-1,1);
    
    % Skip the first index! Since it's the point itself
    % Compute the features regarding the point itself and its neighbors
    for jdx = 2:size(currentNeighbors,2)
        % Normal vector of the current neighbor
        neighborNormalVec = normals(jdx,:);
        % Coordinate of the current neighbor
        neighborCoordinate = coords(jdx,:);
        % Need to select 'source' and 'target' according to the angles
        % between the normal vectors and the line that connects the points
        if dot(currentNormalVec,neighborCoordinate - currentCoordinate) ...
                <= dot(neighborNormalVec,currentCoordinate - neighborCoordinate)
            ps = currentCoordinate;
            pt = neighborCoordinate;
            % Same thing as u = n_{i/j}
            ni = currentNormalVec;
            u = ni;
            nj = neighborNormalVec;
        else
            ps = neighborCoordinate;
            pt = currentCoordinate;
            % Same thing as u = n_{i/j}
            ni = neighborNormalVec;
            u = ni;
            nj = currentNormalVec;
        end
        % Having decided the identities of source and target points we can 
        % compute the 'Darboux frame'
        % Here: (pt = pj) and (ps = pi)
        v = cross(pt - ps,u);
        w = cross(u,v);
        % Now we can ACTUALLY COMPUTE THE FEATURES
        alpha = dot(v,nj);
        phi = dot(u,pt - ps)/norm(pt - ps);
        theta = atan2(dot(w,nj),dot(u,nj));
        % Save the features in the matrix for faster computation of the
        % histograms
        alphas(jdx-1,:) = alpha;
        phis(jdx-1,:)   = phi;
        thetas(jdx-1,:) = theta;
    end
    % Up until now, for point 'idx' all the features regarding itself and
    % its neighbors have been computed.
    % Next thing is to compute the histograms:
    fts(idx,:) = [histcounts(alphas,alphaEdges,'Normalization','probability'),...
        histcounts(phis,phiEdges,'Normalization','probability'),...
        histcounts(thetas,thetaEdges,'Normalization','probability')];
end
fprintf('Finished computing SPFH. Time: %g\n',toc(t1))
% When this for-loop has concluded, each point in the cloud has all its
% features computed (the Simplified Point Feature Histogram (SPFH))

% For computing the FPFH we need to 'then in a second step, for each point 
% we re-determine its k neighbors and use the neighboring SPFH values to 
% weight the final histogram of p (called FPFH):
% FPFH(p) = SPFH(p) + 1/k*\sum_{i=1}^{k}(SPFH(p_{i})/w_{i})
% where the weight w_{i} represents the distance between query point p and 
% a neighbor point p_{i} in a given metric space.

% First, preallocate memory for the FPFH
fpfh = zeros(size(coords,1),3*bins);

% Again, iterate through all the points and compute the actual FPFH
fprintf('Beginning to compute FPFH\n')
for idx = 1:size(coords,1)
    if mod(idx,1000) == 0
        fprintf('Iteration %g of %g\n',idx,size(coords,1))
    end
    % Extract indices of neighbors
    indices = neighbors{idx}(2:end);
    % Extract distances to those neighbors
    currentDistances = dist{idx}(2:end);
    % Extract features of those neighbors
    currentFts = fts(indices,:);
    % k == number of histograms being summed
    k = size(currentFts,1);
    % Compute division (SPFH(p_{i})/w_{i}), sum and division by 'k', which
    % is the second term of the sum
    term = bsxfun(@times,1./(currentDistances),currentFts');
    secondTerm = 1/k*sum(term,2);
    % Sum with the original SPFH and store in 'fpfh'
    fpfh(idx,:) = fts(idx,:) + secondTerm';
    % Normalize (every group of features should add up to 100.0)
    for jdx = bins*(0:2) + 1
        subfts = fpfh(idx,jdx:(jdx+bins - 1));
        fpfh(idx,jdx:(jdx+bins - 1)) = 100*subfts/sum(subfts);
    end
end
fprintf('Finished computing FPFH\n')