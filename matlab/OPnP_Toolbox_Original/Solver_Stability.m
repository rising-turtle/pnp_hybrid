function Solver_Stability(rotation_type, point_config)
%the script to test the numerical stability of various solvers
%including:
%   'Symmetric GB' -- the Grobner basis solver using 2-fold symmetry
%   'Symmetric GB+Polish' - the above solver with one-step DN polishing
%   'Blind GB' -- the Grobner basis solver without using symmetry
%   'DLS RT' -- the resultant solver of DLS
%Input:
%   rotation_type: 'Fully_Random', 'Near_Cayley_Degenerate' or 'Cayley_Degenerate'
%   point_config: 'Ordinary_3D', 'Quasi_Singular' or 'Planar'
%Usage: 
%example: Solver_Stability('Fully_Random','Ordinary_3D')

addpath OPnP;
addpath(genpath('rpnp1.0'));

if nargin < 2
    point_config = 'Ordinary_3D';
end

if nargin < 1
    rotation_type = 'Fully_Random';
end

% the number of independent trials
n_runs = 5000;

err_list = zeros(n_runs,1); %Symmetric GB
err_list2 = zeros(n_runs,1);%Symmetric GB+Polish
err_list3 = zeros(n_runs,1);%Blind GB
err_list4 = zeros(n_runs,1);%DLS RT

for i = 1:n_runs

    % generate random problems
    [var e1 e2 e3 e4 c1 c2 c3 Q q] = Generate_Random_Data_Full(rotation_type,point_config);
    e1 = e1/norm(e1); e2 = e2/norm(e2); e3 = e3/norm(e3); e4 = e4/norm(e4);
    c1 = c1/norm(c1); c2 = c2/norm(c2); c3 = c3/norm(c3); 

    E = [e1;e2;e3;e4];
    C = [c1;c2;c3];

    problem_set(i).groundtruth = var';
    problem_set(i).C = C;
    problem_set(i).E = E;
    problem_set(i).Q = Q;
    problem_set(i).q = q;
end

ids = 1:n_runs;

disp('------ Two-Fold Symmetry GB solver --------');
tic;
for i = ids
    % reorder the monomials
    E = problem_set(i).E;
    
    % solve the equations
    [x y z t] = GB_Solver_3Order_4Variable_Symmetry(E(1,:),E(2,:), E(3,:), E(4,:));
        
    % evaluate errors
    sols = [x; y; z; t];
    
    %normalize solution to unit-norm
    for j = 1:size(sols,2)
        if norm(sols(:,j)) < 1e-6
            continue;
        else
            sols(:,j) = sols(:,j)/norm(sols(:,j));
        end
    end
    
    % evaluate the errors
    [err]     =  evaluate_solutions(sols,problem_set(i).groundtruth);
    err_list(i) = err;
end
toc;

disp('------ Two-Fold Symmetry GB solver with Polish --------');
tic;
for i = ids
    E = problem_set(i).E;
    Q = problem_set(i).Q;
    q = problem_set(i).q;
    
    % solve the equations
    [x y z t] = GB_Solver_3Order_4Variable_Symmetry(E(1,:),E(2,:), E(3,:), E(4,:));
        
    %polish solutions
    for jj = 1:length(x)
        [x(jj) y(jj) z(jj) t(jj)] = PnP_Polish(Q,q,[x(jj) y(jj) z(jj) t(jj)]);
    end
 
    % evaluate errors
    sols = [x; y; z; t];
        
    %normalize solution to unit-norm
    for j = 1:size(sols,2)
        if norm(sols(:,j)) < 1e-6
            continue;
        else
            sols(:,j) = sols(:,j)/norm(sols(:,j));
        end
    end
    
    % evaluate the errors
    [err]     =  evaluate_solutions(sols,problem_set(i).groundtruth);
    err_list2(i) = err;  
    
end
toc;

disp('------- Blind GB solver -------');

tic;
for i = ids
    % reorder the monomials
    E = problem_set(i).E;

    % solve the equations
    [x y z t] = GB_Solver_3Order_4Variable_b_Division(E(1,:),E(2,:), E(3,:), E(4,:));

    % evaluate errors
    sols2 = [x; y; z; t];

    %normalize solution to unit-norm
    for j = 1:size(sols2,2)
        if norm(sols2(:,j)) < 1e-6
            continue;
        else
            sols2(:,j) = sols2(:,j)/norm(sols2(:,j));
        end
    end
    [err2] = evaluate_solutions(sols2,problem_set(i).groundtruth);
    err_list3(i) = err2;
end
toc;


disp('------- Resultant solver in DLS-------');

tic;
for i = ids
    % reorder the monomials
    C = problem_set(i).C;

    % solve the equations
    [y z t] = Resultant_Solver_DLS(C(1,:),C(2,:), C(3,:));

    x = -ones(1,length(y));
    % evaluate errors
    sols3 = [x; y; z; t];

    %normalize solution to unit-norm
    for j = 1:size(sols3,2)
        if norm(sols3(:,j)) < 1e-6
            continue;
        else
            sols3(:,j) = sols3(:,j)/norm(sols3(:,j));
        end
    end

    [err3] = evaluate_solutions(sols3,problem_set(i).groundtruth);
    err_list4(i) = err3;
end
toc;

%% err distribution
%Since the polish step will give extremely accurate solution, thus some
%of the errors might be zero. We clamp such zeros to the minimum non-zero value
% to make the drawing possible. 
temp = err_list2; temp(temp<=0) = [];
minerror = min(temp);

err_list2(err_list2<=0) = minerror;

[hh,bb] = hist(log10(err_list),20);
[hh2,bb2] = hist(log10(err_list2),20);
[hh3,bb3] = hist(log10(err_list3),20);
[hh4,bb4] = hist(log10(err_list4),20);

figure('color','w','position',[100 100 360 320]);
hold all;
box on;

plot(bb,hh,'r-','linewidth',3);hold on;
plot(bb2,hh2,'b-','linewidth',3); hold on;
plot(bb3,hh3,'g-','linewidth',3); hold on;
plot(bb4,hh4,'m-','linewidth',3); hold on;

xlim([-18 0]);
xtick= -18:2:0;
set(gca,'xtick',xtick);

%xlabel('Log_{10} Absolute Unit Quaternion Error','FontSize',11);
%ylabel('Number of Counts','FontSize',11);
legend('Symmetric GB','Symmetric GB + Polish', 'Blind GB','DLS Resultant');
%print('-dpdf',['comp_err_figure_',num2str(stats.neq),'_',num2str(stats.nmon),'_',num2str(n_runs),'.pdf']);