function plot_PT_NOISE()

% method = {'EPNP', 'HYBRID', 'UPNP'}
fid = fopen('PT_NOISE_10.log')
d = cell2mat(textscan(fid, repmat('%f', 1, 13), 'headerLines', 2, 'whitespace', ' ', 'collectoutput', 1));

pt = d(:,1); 
EPNP_t = d(:,2);  EPNP_r = d(:,4); 
HYBRID_t = d(:,6); HYBRID_r = d(:,8); 
UPNP_t = d(:,10); UPNP_r = d(:,12); 
fclose(fid); 

% method = {'LHM', 'EPnP+GN', 'RPnP', 'DLS', 'DLS+++',  'SDP', 'OPnP', 'OPnP+LM'};

fid = fopen('PT_NOISE_10_MAT.txt')
d = cell2mat(textscan(fid, repmat('%f', 1, 33), 'headerLines', 2, 'whitespace', ' ', 'collectoutput', 1));

LHM_t = d(:,2);  LHM_r = d(:,4); 
EPNP_GN_t = d(:,6);  EPNP_GN_r = d(:,8); 
RPNP_t = d(:,10);  RPNP_r = d(:,12); 
DLS_t = d(:,14); DLS_r = d(:,16); 
DLS_plus_t = d(:,18);  DLS_plus_r = d(:,20); 
SDP_t = d(:,22);  SDP_r = d(:,24); 
OPNP_t = d(:,26); OPNP_r = d(:,28); 
OPNP_LM_t = d(:,30);  OPNP_LM_r = d(:,32); 

fclose(fid); 

%% plot translation figure 
close all;
yrange= [0 2];
figure('color', 'w'); 

box('on');
hold('all');

plot(pt, LHM_t, 'marker', 'x', 'color', 'c', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, EPNP_t, 'marker', 'o', 'color', 'm', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, EPNP_GN_t, 'marker', 'd', 'color', 'y', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, RPNP_t, 'marker', '^', 'color', 'r', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, UPNP_t, 'marker', '+', 'color', 'g', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, DLS_t, 'marker', '>', 'color', 'b', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, DLS_plus_t, 'marker', 's', 'color', 'c', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, SDP_t, 'marker', 'v', 'color', 'm', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, OPNP_t, 'marker', '<', 'color', 'r', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, OPNP_LM_t, 'marker', 'p', 'color', 'g', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, HYBRID_t, 'marker', '*', 'color', 'k', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
ylim(yrange);
xlim(pt([1 end]));
set(gca,'xtick',pt);

% title('PnP Comparison','FontSize',12,'FontName','Arial');
xlabel('Number of features with depth measurement (N)','FontSize',11);
ylabel('translational error (%)','FontSize',11);
%legend(p,legendposition);
legend('LHM', 'EPnP', 'DLS', 'DLS+++', 'SDP', 'OPnP', 'OPnP+LM', 'Hybrid-PnP');
% legend('LHM', 'EPnP', 'EPnP+GN', 'RPnP', 'UPnP', 'DLS', 'DLS+++', 'SDP', 'OPnP', 'OPnP+LM', 'Hybrid');
saveas(gcf,'Translation.png');


%% plot rotation figure
figure('color', 'w'); 

box('on');
hold('all');

plot(pt, LHM_r, 'marker', 'x', 'color', 'c', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
plot(pt, EPNP_r, 'marker', 'o', 'color', 'm', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
% plot(pt, EPNP_GN_r, 'marker', 'd', 'color', 'y', 'markerfacecolor', 'n', 'LineWidth',2,'MarkerSize',8);
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
legend('LHM', 'EPnP', 'DLS', 'DLS+++', 'SDP', 'OPnP', 'OPnP+LM', 'Hybrid-PnP');
% legend('LHM', 'EPnP', 'EPnP+GN', 'RPnP', 'UPnP', 'DLS', 'DLS+++', 'SDP', 'OPnP', 'OPnP+LM', 'Hybrid');
saveas(gcf,'Rotation.png');
end

