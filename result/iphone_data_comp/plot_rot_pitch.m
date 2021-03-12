function plot_rot_pitch()

% method = {'EPNP', 'HYBRID', 'UPNP'}
fid = fopen('rot_pitch_15_hybrid.txt')
d = cell2mat(textscan(fid, repmat('%f', 1, 3), 'headerLines', 0, 'whitespace', ' ', 'collectoutput', 1));

pt = d(:,1); 
HYBRID_r = d(:,2);  
fclose(fid); 

% method = {'LHM', 'EPnP+GN', 'RPnP', 'DLS', 'DLS+++',  'SDP', 'OPnP', 'OPnP+LM'};

fid = fopen('rot_pitch_15.txt');
st = 3;
d = cell2mat(textscan(fid, repmat('%f', 1, 17), 'headerLines', 2, 'whitespace', ' ', 'collectoutput', 1));
LHM_r = d(st:end,2); 
EPNP_GN_r = d(st:end,4); 
RPNP_r = d(st:end,6); 
DLS_r = d(st:end,8); 
DLS_plus_r = d(st:end,10); 
SDP_r = d(st:end,12); 
OPNP_r = d(st:end,14); 
OPNP_LM_r = d(st:end,16); 

fclose(fid); 

%% plot translation figure 
close all;
yrange= [0 2];

box('on');
hold('all');

%% plot rotation figure
% figure('color', 'w'); 

box('on');
hold('all');

plot(pt, LHM_r, 'marker', 'x', 'color', 'c', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, EPNP_r, 'marker', 'o', 'color', 'm', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, EPNP_GN_r, 'marker', 'd', 'color', 'm', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, RPNP_r, 'marker', '^', 'color', 'r', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, UPNP_r, 'marker', '+', 'color', 'g', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, DLS_r, 'marker', '>', 'color', 'b', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, DLS_plus_r, 'marker', 's', 'color', 'c', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, SDP_r, 'marker', 'v', 'color', 'm', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, OPNP_r, 'marker', '<', 'color', 'r', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, OPNP_LM_r, 'marker', 'p', 'color', 'g', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, HYBRID_r, 'marker', '*', 'color', 'k', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
ylim(yrange);
xlim(pt([1 end]));
set(gca,'xtick',pt);

% title('PnP Comparison','FontSize',12,'FontName','Arial');
xlabel('Number of features with depth measurement (N)','FontSize',11);
ylabel('rotational error (degree)','FontSize',11);
%legend(p,legendposition);
legend('LHM', 'EPnP+GN', 'DLS', 'DLS+++', 'SDP', 'OPnP', 'OPnP+LM', 'Hybrid-PnP');
% legend('LHM', 'EPnP', 'EPnP+GN', 'RPnP', 'UPnP', 'DLS', 'DLS+++', 'SDP', 'OPnP', 'OPnP+LM', 'Hybrid');
saveas(gcf,'Rotation_rot_pitch_15.png');
end

