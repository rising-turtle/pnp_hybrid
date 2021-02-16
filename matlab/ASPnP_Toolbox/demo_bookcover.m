%addpath 
addpath data
addpath ASPnP

%load data 
load model_bookcover

%inliers-3D points
U = model_bookcover.fpoint(1:3,(model_bookcover.mask>0));

%inliers-2D points
u = model_bookcover.fpoint(6:7,(model_bookcover.mask>0));

%estimate camera pose using ASPnP
[R t] = ASPnP(U,u,model_bookcover.K);

%calculate the projection of vertices
v2d_i = model_bookcover.K*(R*model_bookcover.templateV3D+t*ones(1,4));
v2d_i = v2d_i./repmat(v2d_i(end,:),3,1);
v2d_i = v2d_i(1:2,:);

%draw image
ShowInlierOutlier(model_bookcover.fpoint(4:5,:),model_bookcover.fpoint(6:7,:),....
                  model_bookcover.templateimage,model_bookcover.inputimage,...
                  model_bookcover.mask,model_bookcover.templateV2D,v2d_i);