function [ F draw_upper draw_lower] = FundamentalMatrix_Sampson_Bilinear_SeDuMi_ScaleUP(point1,point2,tol,maxiter)
%estimate the projective fundamental matrix from point correspondences.
%point1'*F*point2 = 0
%minimize the sampson error, first order approximation to the optimal
%critera
%For trilinear term, use iterated bilinear relaxation;
tic
fprintf('\n\n******** Starting Branch&Bound Algorithm For Fundamental Matrix Estimation********\n\n');

if nargin<4 | isempty(maxiter),
    maxiter=2000;
end

if nargin<3 | isempty(tol),
    tol=1-0.95;
end

nbrpoints=size(point1,2); %number of points

if nbrpoints ~= size(point2,2);
    fprintf('Unmatched Point Correspondence');
    exit;
end

if nbrpoints < 8;
    fprintf('At least 8 point correspondences are needed');
    exit;
end
%epsdiff=1e-7;
epsdiff=0.0;

%translate image centroid to origine (and rescale coordinates)  ut
if size(point1,1)<3
    point1 = [point1;ones(1, nbrpoints)];
    point2 = [point2;ones(1, nbrpoints)];
end

[point1, T1] = normalise2dpts(point1);
[point2, T2] = normalise2dpts(point2);

%initial branching region
rect_xL = [ -1;-1;-1; -1;-1;-1;-1; -1; 0];
rect_xU = [ 1 ;1; 1;  1; 1;  1; 1;  1;  1];

%the constant part of relaxation
[b, At,c, At_l,c_l, Q] = Relaxation_ConstantPart(point1,point2);

%an initial solution from linear method 
[F_n8p] = fundmatrix_n8p(point2,point1); %point1'*F*point2 = 0
%[F_EFNS] = Sampson_EFNS (point1,point2);  %point1'*F*point2 = 0
[e_n8p] = EvaluateSampsonError(point1,point2,F_n8p);
%[e_EFNS] = EvaluateSampsonError(point1,point2,F_EFNS);
bestup = e_n8p;
f = F_n8p;

[hh, rect_xL,rect_xU]=Relaxation_LB_loop(point1,point2,rect_xL,rect_xU,b,At,c,At_l,c_l,Q,'N',bestup); %bounding;
[rect_xL rect_xU]

if bestup > sum(hh.res)
    f = hh.f;
    vopt=sum(hh.res);
else
    vopt = bestup;
end

rect_LB=hh.lowerbound;
rect={hh};

draw_upper = [];
draw_lower = [];

iter=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iter <= maxiter, %Branch and Bound-loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [vk,vkindex]=min(rect_LB); %vk: the minimal lower bounding; to branch
   
   draw_upper = [draw_upper vopt];
   draw_lower = [draw_lower vk];
   
   if bestup > vopt
       bestup = vopt;
   end
   
   vdiff=(vopt-vk); %converge or not/absolute difference
   perc=(vopt-vk)/vopt; %relative  difference
   
   disp(['Iter: ',num2str(iter),' Residual: ',num2str(vopt),' Approximation gap: ',num2str(perc*100),'% Regions: ',num2str(length(rect))]);

   if vdiff<epsdiff || perc<tol,
       %'voila'
       break;
   end
   
   %branch on vkindex
   h=rect{vkindex}; %rect: all the current boxes; 

   [slask,pp]=max(rect_xU(1:9,vkindex)-rect_xL(1:9,vkindex)); %largest interval
   
   tmpxL=rect_xL(pp,vkindex);
   tmpxU=rect_xU(pp,vkindex);
   
   %branching strategy....
   
%    bestsol=h.lambda(pp);
%    alfa=0.2; %minimum shortage of interval relative original
%    if (bestsol-tmpxL)/(tmpxU-tmpxL)<alfa,
%        newborder=tmpxL+(tmpxU-tmpxL)*alfa;
%    elseif (tmpxU-bestsol)/(tmpxU-tmpxL)<alfa;
%        newborder=tmpxU-(tmpxU-tmpxL)*alfa;
%    else
%        newborder=bestsol;
%    end
   
   
%   bisect=(tmpyU+tmpyL)/2;
%   bestsol=h.lambda(pp);
%   alfa=0.8;
%   newborder=bestsol*alfa+bisect*(1-alfa);
   
   newborder=(tmpxU+tmpxL)/2; %bisection
   
%       newborder=h.lambda(pp); %best solution


   curr_xL1=rect_xL(:,vkindex);
   curr_xU1=rect_xU(:,vkindex);
   
   curr_xL2=curr_xL1;
   curr_xU2=curr_xU1;
   
   curr_xU1(pp)=newborder;
   curr_xL2(pp)=newborder;
      
   [h1, curr_xL1,curr_xU1]=Relaxation_LB_loop(point1,point2,curr_xL1,curr_xU1,b,At,c,At_l,c_l,Q,'N',bestup);
   [h2, curr_xL2,curr_xU2]=Relaxation_LB_loop(point1,point2,curr_xL2,curr_xU2,b,At,c,At_l,c_l,Q,'N',bestup);
   
   rect_xL=[rect_xL(:,1:vkindex-1),curr_xL1,curr_xL2,rect_xL(:,vkindex+1:end)];
   rect_xU=[rect_xU(:,1:vkindex-1),curr_xU1,curr_xU2,rect_xU(:,vkindex+1:end)];
   
   vopt1=sum(h1.res); %optimal values for the 1st half
   vopt2=sum(h2.res); %iptimal values for the 2nd half
   
   rect={rect{1:vkindex-1},h1,h2,rect{vkindex+1:end}}; %existing boxes
   rect_LB=[rect_LB(1:vkindex-1),h1.lowerbound,h2.lowerbound,rect_LB(vkindex+1:end)]; %lower bounds for all exiting boxes
   
   if vopt1<vopt, %find the better upper bounding
       vopt=vopt1;
       f=h1.f;
   end
   if vopt2<vopt,
       vopt=vopt2;
       f=h2.f;
   end
   
   %screen and remove useless regions
   removeindex=[];
   for ii=1:length(rect);
       if rect{ii}.lowerbound>vopt,
           %remove!
           removeindex(end+1)=ii;
       end
   end
   rect(removeindex)=[];
   rect_xL(:,removeindex)=[];
   rect_xU(:,removeindex)=[];
   rect_LB(removeindex)=[];
       
   iter=iter+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end, %Branch and Bound-loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%denormalize
F = T2'*f*T1;

fprintf('******** Ending (L2,L2) Resectioning Algorithm ********\n\n');
toc
return

function [hh, rect_xL, rect_xU]=Relaxation_LB_loop(point1,point2,rect_xL,rect_xU,b,At,c,At_l,c_l,Q,label,bestup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZATION - over one region with convex envelope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbrpoints=size(point1,2);

%variable order:
% 9 variables in F: f1,f2,....,f9;
% 9 squared terms: f1^2, f2^2, ...,f9^2,
% 12 bilinear terms in denominator: f1f2, f1f3, f2f3,...,f2f5, f2f8,f5f8;
% nbrpoints auxiliary variables 
% 6 trilinear terms in rank constraint;
% 18 dummy bilinear variables

vars = 54 + nbrpoints;

clear K;
K.f = 2;
K.l=0;
K.q = Q;

feasible=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convex relaxation for quadratic terms, socp+linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:9
    %socp part
    At_temp=sparse(zeros(3,vars)); %cone constraints
    c_temp=sparse(zeros(3,1));
    
    At_temp(1,9+i) = 1;
    At_temp(2,9+i) = 1;
    At_temp(3,i) = 2;
    
    c_temp = [1;-1;0];
    
    At = [At; At_temp];
    c = [c; c_temp];
    K.q = [K.q, 3];
    
    %linear part
    At_ltemp=sparse(zeros(1,vars)); %linear inequalities
    c_ltemp=sparse(zeros(1,1));
    
    At_ltemp(i) = (rect_xU(i)+rect_xL(i));
    At_ltemp(9+i) = -1;
    c_ltemp = -rect_xU(i)*rect_xL(i);    
    
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convex relaxation for bilinear terms in denominator, linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VarIndex = [1 2 19;1 3 20; 2 3 21;...
            4 5 22;4 6 23; 5 6 24;...
            1 4 25;1 7 26; 4 7 27;...
            2 5 28;2 8 29; 5 8 30];
for i = 1:12
    At_ltemp=sparse(zeros(4,vars)); %linear inequalities
    c_ltemp=sparse(zeros(4,1));    
    
    %convex envelope
    At_ltemp(1,VarIndex(i,3)) = 1;
    At_ltemp(2,VarIndex(i,3)) = 1;
            
    At_ltemp(1,VarIndex(i,1)) = - rect_xL(VarIndex(i,2));
    At_ltemp(1,VarIndex(i,2)) = - rect_xL(VarIndex(i,1));
    c_ltemp(1) = rect_xL(VarIndex(i,2))*rect_xL(VarIndex(i,1));
        
    At_ltemp(2,VarIndex(i,1)) = - rect_xU(VarIndex(i,2));
    At_ltemp(2,VarIndex(i,2)) = - rect_xU(VarIndex(i,1));
    c_ltemp(2) = rect_xU(VarIndex(i,2))*rect_xU(VarIndex(i,1));
        
    %concave envelope
    At_ltemp(3,VarIndex(i,3)) = -1;
    At_ltemp(4,VarIndex(i,3)) = -1;
        
    At_ltemp(3,VarIndex(i,1)) = rect_xU(VarIndex(i,2));
    At_ltemp(3,VarIndex(i,2)) = rect_xL(VarIndex(i,1));
    c_ltemp(3) = -rect_xU(VarIndex(i,2))*rect_xL(VarIndex(i,1));
        
    At_ltemp(4,VarIndex(i,1)) = rect_xL(VarIndex(i,2));
    At_ltemp(4,VarIndex(i,2)) = rect_xU(VarIndex(i,1));
    c_ltemp(4) = -rect_xL(VarIndex(i,2))*rect_xU(VarIndex(i,1)); 
        
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];
end           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convex relaxation for trilinear terms, (18 auxiliary bilinear terms)linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VarIndex = [1 5 37+nbrpoints;1 9 38+nbrpoints; 5 9 39+nbrpoints;...
            1 6 40+nbrpoints;1 8 41+nbrpoints; 6 8 42+nbrpoints;...
            2 4 43+nbrpoints;2 9 44+nbrpoints; 4 9 45+nbrpoints;...
            2 6 46+nbrpoints;2 7 47+nbrpoints; 6 7 48+nbrpoints;...
            3 4 49+nbrpoints;3 8 50+nbrpoints; 4 8 51+nbrpoints;...
            3 5 52+nbrpoints;3 7 53+nbrpoints; 5 7 54+nbrpoints];
for i = 1:18
    At_ltemp=sparse(zeros(4,vars)); %linear inequalities
    c_ltemp=sparse(zeros(4,1));    
    
    %convex envelope
    At_ltemp(1,VarIndex(i,3)) = 1;
    At_ltemp(2,VarIndex(i,3)) = 1;
            
    At_ltemp(1,VarIndex(i,1)) = - rect_xL(VarIndex(i,2));
    At_ltemp(1,VarIndex(i,2)) = - rect_xL(VarIndex(i,1));
    c_ltemp(1) = rect_xL(VarIndex(i,2))*rect_xL(VarIndex(i,1));
        
    At_ltemp(2,VarIndex(i,1)) = - rect_xU(VarIndex(i,2));
    At_ltemp(2,VarIndex(i,2)) = - rect_xU(VarIndex(i,1));
    c_ltemp(2) = rect_xU(VarIndex(i,2))*rect_xU(VarIndex(i,1));
        
    %concave envelope
    At_ltemp(3,VarIndex(i,3)) = -1;
    At_ltemp(4,VarIndex(i,3)) = -1;
        
    At_ltemp(3,VarIndex(i,1)) = rect_xU(VarIndex(i,2));
    At_ltemp(3,VarIndex(i,2)) = rect_xL(VarIndex(i,1));
    c_ltemp(3) = -rect_xU(VarIndex(i,2))*rect_xL(VarIndex(i,1));
        
    At_ltemp(4,VarIndex(i,1)) = rect_xL(VarIndex(i,2));
    At_ltemp(4,VarIndex(i,2)) = rect_xU(VarIndex(i,1));
    c_ltemp(4) = -rect_xL(VarIndex(i,2))*rect_xU(VarIndex(i,1)); 
        
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];
end           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bounds on the 18 auxiliary bilinear terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:18
    temp_v = [rect_xL(VarIndex(i,1))*rect_xL(VarIndex(i,2))
              rect_xL(VarIndex(i,1))*rect_xU(VarIndex(i,2))
              rect_xU(VarIndex(i,1))*rect_xL(VarIndex(i,2))
              rect_xU(VarIndex(i,1))*rect_xU(VarIndex(i,2))];
    auxiL(i,1) = min(temp_v);
    auxiU(i,1) = max(temp_v);
    
    At_ltemp=sparse(zeros(2,vars)); %linear inequalities
    c_ltemp=sparse(zeros(2,1));
    
    At_ltemp(1,VarIndex(i,3)) = 1;
    At_ltemp(2,VarIndex(i,3)) = -1;
    c_ltemp(1) = -auxiL(i,1);
    c_ltemp(2) = auxiU(i,1);
    
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relaxation for 18 trilinear terms through bilinear relaxation
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
VarIndex = [9 37+nbrpoints 31+nbrpoints;5 38+nbrpoints 31+nbrpoints; 1 39+nbrpoints 31+nbrpoints;...
            8 40+nbrpoints 32+nbrpoints;6 41+nbrpoints 32+nbrpoints; 1 42+nbrpoints 32+nbrpoints;...
            9 43+nbrpoints 33+nbrpoints;4 44+nbrpoints 33+nbrpoints; 2 45+nbrpoints 33+nbrpoints;...
            7 46+nbrpoints 34+nbrpoints;6 47+nbrpoints 34+nbrpoints; 2 48+nbrpoints 34+nbrpoints;...
            8 49+nbrpoints 35+nbrpoints;4 50+nbrpoints 35+nbrpoints; 3 51+nbrpoints 35+nbrpoints;...
            7 52+nbrpoints 36+nbrpoints;5 53+nbrpoints 36+nbrpoints; 3 54+nbrpoints 36+nbrpoints];
for i = 1:18
    At_ltemp=sparse(zeros(4,vars)); %linear inequalities
    c_ltemp=sparse(zeros(4,1));    
    
    %convex envelope
    At_ltemp(1,VarIndex(i,3)) = 1;
    At_ltemp(2,VarIndex(i,3)) = 1;
            
    At_ltemp(1,VarIndex(i,1)) = - auxiL(i);
    At_ltemp(1,VarIndex(i,2)) = - rect_xL(VarIndex(i,1));
    c_ltemp(1) = auxiL(i)*rect_xL(VarIndex(i,1));
        
    At_ltemp(2,VarIndex(i,1)) = - auxiU(i);
    At_ltemp(2,VarIndex(i,2)) = - rect_xU(VarIndex(i,1));
    c_ltemp(2) = auxiU(i)*rect_xU(VarIndex(i,1));
        
    %concave envelope
    At_ltemp(3,VarIndex(i,3)) = -1;
    At_ltemp(4,VarIndex(i,3)) = -1;
        
    At_ltemp(3,VarIndex(i,1)) = auxiU(i);
    At_ltemp(3,VarIndex(i,2)) = rect_xL(VarIndex(i,1));
    c_ltemp(3) = -auxiU(i)*rect_xL(VarIndex(i,1));
        
    At_ltemp(4,VarIndex(i,1)) = auxiL(i);
    At_ltemp(4,VarIndex(i,2)) = rect_xU(VarIndex(i,1));
    c_ltemp(4) = -auxiL(i)*rect_xU(VarIndex(i,1)); 
        
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];
end           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bounds on variables F, linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:9
    At_ltemp=sparse(zeros(2,vars)); %linear inequalities
    c_ltemp=sparse(zeros(2,1));
    
    At_ltemp(1,i) = 1;
    At_ltemp(2,i) = -1;
    c_ltemp(1) = -rect_xL(i);
    c_ltemp(2) = rect_xU(i);
    
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];
end

K.l = size(At_l,1)-K.f;
pars=[];
pars.fid=0;
pars.eps = 1e-8;
% SEDUMI
[x,y,info]=sedumi(-[At_l;At]',b,[c_l;c],K,pars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of OPTIMIZATION - over one region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  feasible==0 || info.dinf==1, %no feasible solution
    res=Inf*ones(1,nbrpoints); %set residuals to infinity
    lowerbound=Inf;
    upperbound = Inf;
    f_initial = [];
else
    %feasible solution
    f_initial = y(1:9);
    f_initial = (1/norm(f_initial))*f_initial;
    f_initial = reshape(f_initial,3,3);
    f_initial = f_initial';
    [U,D,V] = svd(f_initial,0);
    f_initial = U*diag([D(1,1) D(2,2) 0])*V';
    %compute lower bound       
    lowerbound=sum(y([31:30+nbrpoints]))/10000;
    
    %use SNOPT to calculate a upper bound and refine the f
%     f_initial = f_initial';
%     [f, upperbound] = CurrentUpperBound(f_initial(:),point1,point2);
%     f_initial = reshape(f,3,3)';
%     [U,D,V] = svd(f_initial,0);
%     f_initial = U*diag([D(1,1) D(2,2) 0])*V';
    
    upperbound = 0;
    for i = 1:nbrpoints
        Mat1 = f_initial*point2(:,i);
        Mat2 = f_initial'*point1(:,i);
        upperbound = upperbound + (point1(:,i)'*f_initial*point2(:,i))^2/(Mat1(1)^2+Mat1(2)^2+Mat2(1)^2+Mat2(2)^2);
    end
    
%     %Refine upper bound through fmincon
%     [f, upperbound] = RefineUpperBoundFmincon(f_initial(:),point1,point2);
%     f_initial = reshape(f,3,3)';
%     [U,D,V] = svd(f_initial,0);
%     f_initial = U*diag([D(1,1) D(2,2) 0])*V';
%     
%     upperbound = 0;
%     for i = 1:nbrpoints
%         Mat1 = f_initial*point2(:,i);
%         Mat2 = f_initial'*point1(:,i);
%         upperbound = upperbound + (point1(:,i)'*f_initial*point2(:,i))^2/(Mat1(1)^2+Mat1(2)^2+Mat2(1)^2+Mat2(2)^2);
%     end

%bounds contraction
K.l = K.l + 1;
if label == 'Y'
    %refine solution
    if bestup < upperbound
        upperbound = bestup;
    end    
   
    for i = 1:9
        b_temp = zeros(vars,1);
        b_temp(i) = 1;
        [x,y,info]=sedumi(-[At_l;b';At]',-b_temp,[c_l;upperbound*10000;c],K,pars);
        rect_xL(i) = y(i);
        
        [x,y,info]=sedumi(-[At_l;b';At]',b_temp,[c_l;upperbound*10000;c],K,pars);
        rect_xU(i) = y(i);
    end
end     
 
end
hh.f = f_initial;
hh.res=upperbound; %residuals
hh.lowerbound=lowerbound; %lowerbounds for each fractional term

%the constant part of the relaxation
function [b, At,c,At_l,c_l,Q] = Relaxation_ConstantPart(point1,point2)
nbrpoints=size(point1,2);

%variable order:
% 9 variables in F: f1,f2,....,f9;
% 9 squared terms: f1^2, f2^2, ...,f9^2,
% 12 bilinear terms in denominator: f1f2, f1f3, f2f3,...,f2f5, f2f8,f5f8;
% nbrpoints auxiliary variables 
% 6 trilinear terms in rank constraint;
% 18 dummy bilinear variables

vars = 54 + nbrpoints;

%sedumi matrices
At_l=sparse(zeros(0,vars)); %linear inequalities
c_l=sparse(zeros(0,1));
At=sparse(zeros(0,vars)); %cone constraints
c=sparse(zeros(0,1));
clear K;
K.l=0;
K.q=[];

%objective function
b=sparse(zeros(vars,1));
b(31:30+nbrpoints)=-1; %minimize sum{delta(i)}

At_ltemp=sparse(zeros(1,vars)); %linear equalities
c_ltemp=sparse(zeros(1,1));
%scale ambiguity    
At_ltemp(10:18) = ones(1,9);
c_ltemp = -1;

At_l = [At_l; At_ltemp];
c_l = [c_l; c_ltemp]; 

At_ltemp=sparse(zeros(1,vars)); %linear equalities
c_ltemp=sparse(zeros(1,1));
%rank constraint  
At_ltemp(30+nbrpoints+1:30+nbrpoints+6) = [1,-1,-1,1,1,-1];
c_ltemp = 0;

At_l = [At_l; At_ltemp];
c_l = [c_l; c_ltemp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% residuals 1,2,3,...,nbrpoints, socp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nbrpoints
    At_temp=sparse(zeros(3,vars)); %cone constraints
    c_temp=sparse(zeros(3,1));
    
    At_temp(1,10:17) = [point2(1,i)^2+point1(1,i)^2, point2(2,i)^2+point1(1,i)^2, point2(3,i)^2, ... %%%error
                        point2(1,i)^2+point1(2,i)^2, point2(2,i)^2+point1(2,i)^2, point2(3,i)^2, ...
                        point1(3,i)^2, point1(3,i)^2];
    At_temp(1,19:21) = [2*point2(1,i)*point2(2,i),2*point2(1,i)*point2(3,i),2*point2(2,i)*point2(3,i) ];
    At_temp(1,22:24) = At_temp(1,19:21);
    At_temp(1,25:27) = [2*point1(1,i)*point1(2,i),2*point1(1,i)*point1(3,i),2*point1(2,i)*point1(3,i) ];
    At_temp(1,28:30) = At_temp(1,25:27);
    At_temp(1,30+i) = 1;
    
    At_temp(2,:) = At_temp(1,:);
    At_temp(2,30+i) = -1;
    
    At_temp(3,1:9) = 100*2*point2(:,i)'*[point1(1,i)*eye(3) point1(2,i)*eye(3) point1(3,i)*eye(3)]; %%%error
    
    At = [At; At_temp];
    c = [c; c_temp];
    K.q = [K.q 3];
    
    %the squared denominator should be positive
    At_ltemp=sparse(zeros(1,vars)); %linear inequalities
    c_ltemp=sparse(zeros(1,1));
    
    At_ltemp = At_temp(1,:);
    At_ltemp(30+i) = 0;
    
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];      
end

%the auxiliary variables should be positive
for i = 1:nbrpoints
    At_ltemp=sparse(zeros(1,vars)); %linear inequalities
    c_ltemp=sparse(zeros(1,1));
    
    At_ltemp(30+i) = 1;
        
    At_l = [At_l; At_ltemp];
    c_l = [c_l; c_ltemp];      
end
    
Q = K.q;

function [f upperbound] = CurrentUpperBound(f_initial,point1,point2)
xlow = [-1;-1;-1;-1;-1;-1;-1;-1;0];
xupp = [1;1;1;1;1;1;1;1;1];
Flow = [-inf;1;0];
Fupp = [inf; 1;0];
assignin('base','point1',point1);
assignin('base','point2',point2);
assignin('base','delta1',zeros(size(point1)));
assignin('base','delta2',zeros(size(point2)));
[f,F,inform,xmul,Fmul] = snopt(f_initial,xlow,xupp,Flow,Fupp,'ValueGradient');
upperbound = F(1);

function [f, upperbound] = RefineUpperBoundFmincon(f_initial,point1,point2);
xlow = [-1;-1;-1;-1;-1;-1;0.6;-1;-1];
xupp = [1;1;1;1;1;1;1;1;1];
options = optimset('GradObj','on','Display','off','GradConstr','on');
[x,fval,exitflag] = fmincon(@FminconObj,f_initial,[],[],[],[],xlow,xupp,@FminconConstraint,options);
f = x;
upperbound = fval;










