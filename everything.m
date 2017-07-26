% [pt1, feat1] = read_features('../../dataset/pairwise_noise_xyz_level_02_01_rot_05/features_0000.bin');
% [pt2, feat2] = read_features('../../dataset/pairwise_noise_xyz_level_02_01_rot_05/features_0001.bin');
%% Define parameters
farDistThr          = 1.8;
viewPoint           = [0,0,0];

%% Load point clouds to be merged
normals     = cell(1,4);
coords      = cell(1,4);
ptClouds    = cell(1,4);

figure(1)
for idx = 1:size(ptClouds,2)
    % Load from .txt
    tempData = load(sprintf('mk2cap%g.txt',idx));
    % Remove Inf's, NaN's and points that are far away
    tempData = cleanData(tempData,farDistThr);
    % Take points that tend to red color
    indices = (tempData(:,4) > 90) & (tempData(:,5) < 80) & ...
        (tempData(:,6) < 80);
    tempData = tempData(indices,:);
    % Convert to pointCloud object, denoise, and store it in cell array
    ptClouds{idx} = pcdenoise(pointCloud(double(tempData(:,1:3)),'Color',...
        uint8(tempData(:,4:6))),'Threshold',1e-3);
    % Extract coordinates
    coords{idx} = ptClouds{idx}.Location;
    % Compute normal vectors at each point
    normals{idx} = pcnormals(ptClouds{idx});

    % Compute orientations consistently according to the camera's viewpoint
    where = dot(bsxfun(@minus,viewPoint,coords{idx})', normals{idx}') < 0;
    normals{idx}(where',:) = -1*normals{idx}(where',:);
    % Display de point clouds
    subplot(3,4,idx), pcshow(ptClouds{idx}),...
        title(sprintf('Point cloud %g',idx))
end

%% Compute features
% Load coordinates and normals, since the server doesn't have ptCloud objs
% load('coordsAndNormals.mat')
% fprintf('Loaded coordinates and normal vectors\n')
% Define number of bins per feature (remember: 3 features!)
bins = 11;
% Define radius
radius = .04;% 20 cm??
% Cell array to store the features of the various point clouds
features = cell(1,4);
tt = tic;
fprintf('Starting to compute features\n')
for idx = 1:size(coords,2)
    fprintf('\nFeatures of point cloud number %g\n',idx)
    features{idx} = computeFPFHFeatures(coords{idx},normals{idx},radius,...
        bins);
end
fprintf('Finished computing features. Time: %g\n',toc(tt))

save('myfeats.mat','features','coords','ptClouds')
% Clear variables except the two that are needed
clearvars -except coords features ptClouds

%% Perform fast global registration
cd(['/Users/j1k1000o/Documents/MATLAB/Kinect/FastGlobalRegistration/',...
    'FastGlobalRegistration/source/Matlab'])
transformation = cell(1,4);
transformedCoords = cell(1,4);
transformedCoords{1} = coords{1};
% I'll just assume that everything is to be computed from the viewpoint of
% the first point cloud
sourceLocation = coords{1};
sourceFeatures = features{1};
for idx = 1:4
    % Compute transformation
    fprintf('Computing transformation %g\n',idx)
    trans = tic;
    transformation{idx} = fast_global_registration(...
        sourceLocation,sourceFeatures,...
        coords{idx},features{idx});
    fprintf('Finished computing transformation %g. Time: %g\n',...
        idx,toc(trans))
    disp(transformation{idx})
    % Apply transformation
    fprintf('Applying transformation %g\n',idx)
    transformedCoords{idx} = (transformation{idx}(1 : 3, 1 : 3) * ...
        (coords{idx})' + repmat(transformation{idx}(1 : 3, 4), 1, ...
        size(coords{idx}, 1)))';
    fprintf('Transformation %g was applied\n',idx)
end

%% Generate pointCloud objects for easier display
newPtClouds = cell(1,4);
figure(1)
for idx = 1:4
    newPtClouds{idx} = pointCloud(transformedCoords{idx},'Color',...
        ptClouds{idx}.Color);
    subplot(3,4,idx+4), pcshow(newPtClouds{idx}),...
        title(sprintf('New point cloud %g',idx))
    subplot(3,4,9:12), hold on, pcshow(newPtClouds{idx})
end

