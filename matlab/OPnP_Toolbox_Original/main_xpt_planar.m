clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('gOp'));
addpath dls_pnp_matlab;
addpath OPnP;

% experimental parameters
nls= 0.5:0.5:5;
npt= 10;
num= 400;

% compared methods
A= zeros(size(nls));
B= zeros(num,1);
name= {'SP+LHM', 'EPnP', 'RPnP', 'DLS', 'DLS+++',  'SDP'      'OPnP', 'OPnP+LM'};
f= {    @SP, @EPnP,  @RPnP, @dls_pnp_all, @robust_dls_pnp, @gOp_interface,@OPnP,@PnP_Reproj_NLS_Matlab};
marker= { 'x', 'o', 'd', '^', '+', '>', 's', 'none'};
color= {'r','g','b','m','m', 'y','k','c'};
markerfacecolor=  {'r','g','n','m','n', 'n','n', 'n'};
linestyle= {'-','-','-','-',':','-','-','-.'};

method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor, 'linestyle', linestyle);

% experiments
for i= 1:length(nls)
    
    nl= nls(i);
    fprintf('nl = %.1f: ',nl);
    
    for k= 1:length(method_list)
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    
    index_fail = [];
    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;
        f= 800;
        
        % generate 3d coordinates in camera space
        XXw= [xrand(2,npt,[-2 2]); zeros(1,npt)];
		R= rodrigues(randn(3,1));
        t= [rand-0.5;rand-0.5;rand*8+4];
        Xc= R*XXw+repmat(t,1,npt);
        
        % projection
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;

         % pose estimation
         for k= 1:length(method_list)
             if strcmp(method_list(k).name, 'OPnP+LM')
                 for jj = 1:size(t1,2)
                     [R1(:,:,jj),t1(:,jj)]= method_list(k).f(XXw,xxn/f,R1(:,:,jj),t1(:,jj));
                 end
             else
                try
                    [R1,t1]= method_list(k).f(XXw,xxn/f);
                catch
                    fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                    index_fail = [index_fail, j];
                    break;
                end
             end
            
            %no solution
            if size(t1,2) < 1
                fprintf(['   The solver - ',method_list(k).name,' - returns no solution! \n']);
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
yrange= [0 5];

i= 0; w= 350; h= 350;

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'mean_r','Mean Rotation Error',...
    'Gaussian Image Noise (pixels)','Rotation Error (degrees)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_r','Median Rotation Error',...
    'Gaussian Image Noise (pixels)','Rotation Error (degrees)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'mean_t','Mean Translation Error',...
    'Gaussian Image Noise (pixels)','Translation Error (%)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_t','Median Translation Error',...
    'Gaussian Image Noise (pixels)','Translation Error (%)',2);

rmpath(genpath('rpnp1.0'));
rmpath(genpath('gOp'));
rmpath dls_pnp_matlab;
rmpath OPnP;
