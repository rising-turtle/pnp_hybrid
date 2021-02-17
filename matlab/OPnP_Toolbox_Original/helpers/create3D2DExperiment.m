function [XXw, xxn, t12, R12 ] = create3D2DExperiment( pt_number, noise, cam, outlier_fraction )
% function [v1, v2, d1, t12, R12 ] = create3D2DExperiment( pt_number, noise, cam, outlier_fraction )

%%
%  output 
% v1, v2: bearingsVector 
% d1: depth vector 
% t12, R12: randomly generated transformation 

%% generate random view-points
max_parallax = 2.0;
max_rotation = 0.5;

position1 = zeros(3,1);
rotation1 = eye(3);

position2 = max_parallax * 2.0 * (rand(3,1) - repmat(0.5,3,1));
rotation2 = generateBoundedR(max_rotation);

%% Generate random point-cloud

minDepth = 2.0;
maxDepth = 8.0;

% normalizedPoints = 2.0*(rand(3,pt_number)-repmat(0.5,3,pt_number));
% norms = sqrt(sum(normalizedPoints.*normalizedPoints));
% directions = normalizedPoints./repmat(norms,3,1);
% points = (maxDepth-minDepth) * normalizedPoints + minDepth * directions;

% points = generateRandomPoint(cam, minDepth, maxDepth, pt_number);

% focal_length = 800.0;
focal_length = cam.F; 

v1 = zeros(3,pt_number);
v2 = zeros(3,pt_number);
d1 = zeros(1,pt_number); 

% for i=1:pt_number
i = 1;
while i<= pt_number    
    pt  =  generateRandomPoint(cam, minDepth, maxDepth);
    body_point1 = rotation1' * (pt-position1);
    body_point2 = rotation2' * (pt-position2);
    
    if body_point2(3) <= 0.
        continue;
    end
    
    d_std = poly_curve(body_point1(3));
    % d = body_point1(3) + 2.0*(rand(1) - 0.5)*d_std; 
    
    d_std = 0;
    noise_d = d_std * randn(1);
    d = noise_d + body_point1(3); 
    d1(i) = d; 
    % we actually omit the can rotation here by unrotating the bearing
    % vectors already
    bearingVector1 = body_point1 ;
    bearingVector2 = body_point2 ;
    bearingVector1_norm = norm(bearingVector1);
    bearingVector2_norm = norm(bearingVector2);
    bearingVector1 = bearingVector1/bearingVector1_norm;
    bearingVector2 = bearingVector2/bearingVector2_norm;
    
    % add noise to the bearing vectors here
    bearingVector1_noisy = addNoise(bearingVector1,focal_length,noise);
    bearingVector2_noisy = addNoise(bearingVector2,focal_length,noise);
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector1_norm = norm(bearingVector1_noisy);
    bearingVector2_norm = norm(bearingVector2_noisy);
    
    v1(:,i) = [bearingVector1_noisy./bearingVector1_norm];
    v2(:,i) = [bearingVector2_noisy./bearingVector2_norm];
    i = i + 1;
end

%% Add outliers
number_outliers = floor(outlier_fraction*pt_number);

if number_outliers > 0
for i=1:number_outliers
    
    %generate random point
    normalizedPoint = 2.0*(rand(3,1)-repmat(0.5,3,1));
    norm1 = sqrt(sum(normalizedPoint.*normalizedPoint));
    direction = normalizedPoint./norm1;
    point = (maxDepth-minDepth) * normalizedPoint + minDepth * direction;
    
    body_point2 = rotation2' * (point-position2);
    
    % store the point (no need to add noise)
    bearingVector2 = body_point2;
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector2_norm = norm(bearingVector2);
    
    v2(:,i) = [bearingVector2./bearingVector2_norm];
end
end

%% compute relative translation and rotation

R12 = rotation1' * rotation2;
t12 = rotation1' * (position2 - position1);

R21 = R12';
t21 = -R21 * t12; 
%% Cut the cam offset in the single camera case (e.g. central)

% if cam_number == 1
%    v1 = v1(1:3,:);
%    v2 = v2(1:3,:);
% end

XXw = v1(1:3,:)./v1(3,:);
XXw = XXw(1:3,:).*d1;
xxn = v2(1:2,:)./v2(3,:); 

end

function pt = generateRandomPoint(cam, min_d, max_d)
    pt = zeros(3,1); 
    u = rand(1)*2*cam.CX; 
    v = rand(1)*2*cam.CY; 
    d = rand(1)*(max_d- min_d) + min_d;
    pt(1) = d*(u-cam.CX)/cam.F; 
    pt(2) = d*(v-cam.CY)/cam.F; 
    pt(3) = d; 
end

function d_std = poly_curve(d)
    d_std = (d^2*0.00155816 - 0.00362021*d + 0.00452812);
end

