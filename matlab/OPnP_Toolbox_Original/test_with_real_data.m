%
% Feb. 25 2021, He Zhang, fuyinzh@gmail.com   
%   
%    run PnP methods using iphone's data (tracked features in 7 consecutive images) 
% 

clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('gOp'));
addpath dls_pnp_matlab;
addpath OPnP;
addpath helpers;
addpath SeDuMi_1_3;

%% experimental parameters
% nl= 1; %2;
npts= 4:2:30; %4:1:15;
num= 10; %

% compared methods
A= zeros(size(npts));
B= zeros(num,1);
name= {'LHM', 'EPnP+GN', 'RPnP', 'DLS', 'DLS+++',  'SDP', 'OPnP', 'OPnP+LM'};
f= {    @LHM, @EPnP_GN,  @RPnP, @dls_pnp_all, @robust_dls_pnp, @gOp_interface,@OPnP,@PnP_Reproj_NLS_Matlab};
marker= { 'x', 'o', 'd', '^', '+', '>', 's', 'none'};
color= {'r','g','b','m','m', 'y','k','c'};
markerfacecolor=  {'r','g','n','m','n', 'n','n', 'n'};
linestyle= {'-','-','-','-',':','-','-','-.'};

method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor, 'linestyle', linestyle);

%% data path 
folder = '../../result/feats_tracker_results'; 
move_name = 'rot_pitch';
gt = [0, 3.13, 6.19, 9.25, 12.32, 15.32, 18.44]; % ground truth for rot_pitch, rotate around Z axis 

%% 
index = 6; 
r_angle = (gt(index+1) - gt(1))*pi/180.; 

t = [0, 0, 0]'; 
cs = cos(r_angle); 
ss = sin(r_angle); 

% pitch rotate around Y axis 
R = [cs, 0, ss; 
     0, 1., 0; 
     -ss ,0 , cs]; 

% experiments
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d: ',npt);
    
     for k= 1:length(method_list)
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    
    index_fail = [];
    
    for j= 1:num
        
        %% get 
        filename = join([folder,'/', move_name,'_', string(j)], "");
        [XXw, xxn] = get_data(filename); 
        xxn = xxn(:,(2*(index-1)+1):(2*index)); 
        
        %% get the number of points 
        [XXw, xxn] = get_N(XXw, xxn, npt); 
        XXw = XXw'; 
        xxn = xxn';
        
         % pose estimation
        for k= 1:length(method_list)
             if strcmp(method_list(k).name, 'OPnP+LM')
                 for jj = 1:size(t1,2)
                     [R1(:,:,jj),t1(:,jj)]= method_list(k).f(XXw,xxn,R1(:,:,jj),t1(:,jj));
                 end
             else
                try
                    [R1,t1]= method_list(k).f(XXw,xxn);
                catch
                    fprintf(['    The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                    index_fail = [index_fail, j];
                    break;
                end
             end
            
            %no solution
            if size(t1,2) < 1
                fprintf(['    The solver - ',method_list(k).name,' - returns no solution! \n']);
                index_fail = [index_fail, j];
                break;
            end
            %choose the solution with smallest error 
            error = inf;
            for jjj = 1:size(R1,3)
                tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                if sum(tempy) < error
                    y = tempy;
                    error = sum(tempy);
                end
            end
            
            if strcmp(method_list(k).name,'OPnP') && error > 30
                txt = 1;
            end
            
            method_list(k).r(j)= y(1);
            method_list(k).t(j)= y(2);
        end

        showpercent(j,num);
    end
    fprintf('\n');
    
    % save result
    for k= 1:length(method_list)
        method_list(k).r(index_fail) = [];
        method_list(k).t(index_fail) = [];
        
        method_list(k).mean_r(i)= mean(method_list(k).r);
        method_list(k).mean_t(i)= mean(method_list(k).t);
        method_list(k).med_r(i)= median(method_list(k).r);
        method_list(k).med_t(i)= median(method_list(k).t);
        method_list(k).std_r(i)= std(method_list(k).r);
        method_list(k).std_t(i)= std(method_list(k).t);
    end
end

close all;
yrange= [0 2];

i= 0; w= 300; h= 300;

output_file = join(['./test_results/', move_name, '_', string(index*3),'.txt'],"");
save_results_R(output_file, npts, method_list); 

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'mean_r','Mean Rotation Error',...
    'Number of Points','Rotation Error (degrees)');

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_r','Median Rotation Error',...
    'Number of Points','Rotation Error (degrees)');


function [XXw, xxn] = get_N(XXw, xxn, npt)
    N = size(xxn,1); 
    seq = randperm(N, npt); 
    XXw = XXw(seq, :);
    xxn = xxn(seq, :); 
end
