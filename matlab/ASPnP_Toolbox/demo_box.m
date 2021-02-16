%addpath 
addpath data
addpath ASPnP

%load data 
load model_box

%inliers-3D points
U = model_box.fpoint(1:3,(model_box.mask>0));

%inliers-2D points
u = model_box.fpoint(6:7,(model_box.mask>0));

%estimate camera pose using ASPnP
[R t] = ASPnP(U,u,model_box.K);

%calculate the projection of vertices
v2d_i = model_box.K*(R*model_box.templateV3D+t*ones(1,8));
v2d_i = v2d_i./repmat(v2d_i(end,:),3,1);
v2d_i = v2d_i(1:2,:);

%draw image
ShowInlierOutlier(model_box.fpoint(4:5,:),model_box.fpoint(6:7,:),....
                  model_box.templateimage,model_box.inputimage,...
                  model_box.mask,model_box.templateV2D,v2d_i);