clear; 

t = load('output_cnt_2d_30_no_opt.log');
x = t(:,1);
epnp_no_ba = t(:,5); 
hybridpnp_no_ba = t(:,3);
t = load('output_cnt_2d_30_opt.log'); 
hybridpnp = t(:,3);
epnp_ba = t(:,5); 

plot(x, hybridpnp_no_ba, 'm+-');
hold on; 
plot(x, hybridpnp, 'rs-'); 
hold on; 
plot(x, epnp_no_ba, 'gd-');
hold on; 
plot(x, epnp_ba, 'bx-');
grid on; 
