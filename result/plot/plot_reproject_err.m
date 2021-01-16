clear; 

t = load('output_cnt_2d_30_no_opt.log');
x = t(:,1);
epnp_no_ba = t(:,11); 
hybridpnp = t(:,10);

plot(x, hybridpnp, 'rs-'); 
hold on; 
plot(x, epnp_no_ba, 'gd-');
hold on; 
% plot(x, epnp_ba, 'bx-');
% grid on; 