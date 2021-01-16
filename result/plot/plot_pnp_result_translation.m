clear; 

t = load('output_cnt_2d_30_no_opt.log');
x = t(:,1);
epnp_no_ba = t(:,4); 
hybridpnp_no_ba = t(:,2);
t = load('output_cnt_2d_30_opt.log'); 
hybridpnp = t(:,2);
epnp_ba = t(:,4); 

plot(x, hybridpnp_no_ba, 'm+-');
hold on; 
plot(x, hybridpnp, 'rs-'); 
hold on; 
plot(x, epnp_no_ba, 'gd-');
hold on; 
plot(x, epnp_ba, 'bx-');
grid on; 