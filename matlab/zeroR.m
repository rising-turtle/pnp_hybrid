clear
close all
format short g
% clc
addpath('/Users/jin/Q_Mac/mexopencv');

testNum = 100; %150;
stdIdx = 3;
secNum = 3;
for dnosie = 1% 0.0001:0.001:0.05
    dnosie
    results = zeros(testNum,8);
    results_pnp = zeros(testNum,8);
    gt = [0,0,0];
    pixelNoise = zeros(testNum,2);
    for testId = 1:testNum
        testId;
        data_path = "/Volumes/Mac_bkp/z/rot_yaw_"; % % z/rot_yaw_ % x/rot_roll_  % y/rot_pitch_
        I1 = imread(data_path+num2str(stdIdx)+"/color/"+num2str(350+testId)+".png");
        I1 = rgb2gray(I1);
%         I1 = imresize(I1,[480,640]);

        I_dpt = imread(data_path+num2str(stdIdx)+"/depth/"+num2str(350+testId)+".png");


        tempCell = cv.goodFeaturesToTrack(I1, 'MaxCorners', 800, 'QualityLevel', 0.01, 'MinDistance', 20);
        prevPts = double(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';


        imgId = 250+testId;
        slct = 3;
        I_cur = imread(data_path+num2str(slct)+"/color/"+num2str(imgId)+".png");
        I_cur = rgb2gray(I_cur);
        I_dpt2 = imread(data_path+num2str(slct)+"/depth/"+num2str(imgId)+".png");

        
%         I_cur = imresize(I_cur,[480,640]);

        tempCell = cv.calcOpticalFlowPyrLK(I1, I_cur, prevPts);
        nextPts = double(reshape(cell2mat(tempCell)',[2,length(tempCell)]))';
        % C2 = double(reshape(tempMat',[2,length(tempMat)/2]))';
        idx = find(nextPts(:,1)>0&nextPts(:,1)<1920&nextPts(:,2)>0&nextPts(:,2)<1440);
        nextPts = nextPts(idx,:);
        prevPts = prevPts(idx,:);

        [F, mask] = cv.findFundamentalMat(prevPts, nextPts, 'Method','Ransac','RansacReprojThreshold',2);
        nextPts = nextPts(mask==1,:);
        prevPts = prevPts(mask==1,:);
        pixelNoise(testId,:) = mean(abs(nextPts - prevPts));
        
        
    %     figure,showMatchedFeatures(I1,I_cur,prevPts,nextPts,'montage','PlotOptions',{'r*','g*','y-'});
    %     break
        dd = importdata(data_path+num2str(slct)+"/Frames.txt");
        dd = dd;
        fx = dd(imgId,3);
        fy = dd(imgId,4);
        cx = dd(imgId,5);
        cy = dd(imgId,6);
        K2 = [fx,0,cx;0,fy,cy;0,0,1];

        csvwrite("intrinsic.csv",[fx,fy,cx,cy]);
        dd = importdata(data_path+num2str(stdIdx)+"/Frames.txt");
        dd = dd;
        fx = dd(imgId,3);
        fy = dd(imgId,4);
        cx = dd(imgId,5);
        cy = dd(imgId,6);
        K1 = [fx,0,cx;0,fy,cy;0,0,1];


        mvImg = [];
        fixObj = [];
        mvObj = [];
        
        cswriteTempMv = [];
        cswriteTempFx = [];
        psdepth = 0.2+4*rand(length(prevPts),1)+dnosie*(wgn(length(prevPts),1,0));
        for j = 1:length(prevPts)
            ui = double(prevPts(j,1));
            vi = double(prevPts(j,2));
            xi = round(ui/7.5);
            yi = round(vi/7.5);
            xi = max(xi,1);
            yi = max(yi,1);
            dist_xi = xi - ui/7.5;
            dist_yi = yi - vi/7.5;
            d1 = double(I_dpt(yi, xi)) * 0.001;
%             d1 = psdepth(j);
            
            uj = double(nextPts(j,1));
            vj = double(nextPts(j,2));
            xj = round(uj/7.5);
            yj = round(vj/7.5);
            xj = max(xj,1);
            yj = max(yj,1);
            dist_xj = xj - uj/7.5;
            dist_yj = yj - vj/7.5;
            d2 = double(I_dpt2(yj, xj)) * 0.001;
            
            closeDepth = 0;


            dis_thred = 0.4;
            if (abs(dist_xi)<dis_thred && abs(dist_yi)<dis_thred && abs(dist_xj)<dis_thred && abs(dist_yj)<dis_thred)
                closeDepth = 1;
            end


            if(d1>= 0.2 && d1 <= 4.5 && d2>= 0.2 && d2 <= 4.5 && closeDepth>0) 
%                 d1
                pt_x = d1*(ui-cx)/fx;
                pt_y = d1*(vi-cy)/fy;
                ptfx = [pt_x,pt_y,d1];
                fixObj = [ptfx;fixObj];
                
                pt_x = d2*(uj-cx)/fx;
                pt_y = d2*(vj-cy)/fy;
                ptmv = [pt_x,pt_y,d2];
                mvObj = [ptmv;mvObj];
                mvImg = [double(nextPts(j,:));mvImg];
                cswriteTempFx = [[prevPts(j,1),prevPts(j,2),d1];cswriteTempFx];
                cswriteTempMv = [[nextPts(j,1),nextPts(j,2)];cswriteTempMv];
            end


        end
        csvwrite("fixUVd.csv",mvObj);
        csvwrite("mvUVd.csv",cswriteTempMv);
        tvec = rand(3,1);
        rvec = rand(3,1);
%         [rvec, tvec] = cv.solvePnP(fixObj, mvImg, K2, 'Rvec',rvec,  'Tvec',tvec, 'Method', 'Iterative');
%         [rotm,tvec] = Ransac3d3d(fixObj,mvObj);
        [rvec, tvec, success, inliers] = cv.solvePnPRansac(fixObj(1:100,:), mvImg(1:100,:), K2, 'Rvec',rvec,  'Tvec',tvec,'Method','EPnP');%P3P


        rotm = cv.Rodrigues(rvec)
        1000*tvec'
        rotm2eul(rotm)*180/pi
        break
        results_pnp(testId,1:3) = 1000*tvec;%rotm2eul(rotm)*180/pi;
        axang = rotm2axang(rotm);
        results_pnp(testId,4) = 1000*norm(tvec);%axang(4)*180/pi;
        results_pnp(testId,5:7) = rotm2eul(rotm)*180/pi;
        results_pnp(testId,8) = axang(4)*180/pi;

%         break
    end
%     break
    % disp('mean') 
    memm = mean(abs(results_pnp))
    % disp('std')
    stdmm = std(results_pnp)
    disp('========================')
end