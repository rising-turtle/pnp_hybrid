clear; 
% num_3d_points (e)-epnp (h)-hybrid mt-mean_t st-std_t mr-mean_r sr-std_r
% num_3d_points e_mt e_st e_mr e_sr h_mt h_st h_mr h_sr 

data = load('rot_pitch_1_2.log');

cnt_3d = data(:,1);

e_mt = data(:,2);
e_st = data(:,3); 
e_mr = data(:,4);
e_sr = data(:,5); 

h_mt = data(:,6); 
h_st = data(:,7); 
h_mr = data(:,8); 
h_sr = data(:,9); 

%% draw results 
figure;
% draw translation 
subplot(1,2,1); 
% plot(cnt_3d, hybridpnp_no_ba, 'm+-');
errorbar(cnt_3d, e_mt, e_st, '-rs', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
hold on; 
% plot(x, hybridpnp, 'rs-'); 
errorbar(cnt_3d, h_mt, h_st, '-.bo', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');
title('Translation error against number of 3D points');
xlabel('Number of 3D points');
ylabel('Translation error norm [m]');
   
% plot(x, epnp_no_ba, 'gd-'); 
% plot(x, epnp_ba, 'bx-');
% grid on; 

% draw rotation 
subplot(1,2,2);
errorbar(cnt_3d, e_mr, e_sr, '-rs', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
hold on; 
errorbar(cnt_3d, h_mr, h_sr, '-.o', 'MarkerSize', 5, ...
        'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'blue');

title('Rotation error against number of 3D points');
xlabel('Number of 3D points');
ylabel('Rotation error [rad]');