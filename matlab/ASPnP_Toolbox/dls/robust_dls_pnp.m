function [C_est, t_est, cost, flag] = robust_dls_pnp(p, z)


R = cat(3, rotx(pi/2), roty(pi/2), rotz(pi/2));
t = mean(p,2);

cost = inf;
C_est = [];
t_est = [];
flag = [];
for i = 1:3
    % Make a random rotation
    pp = R(:,:,i) * (p - repmat(t, 1, size(p,2)));
    
    [C_est_i, t_est_i, cost_i, flag_i] = dls_pnp(pp, z);
    
    for j = 1:length(cost_i)
        t_est_i(:,j) = t_est_i(:,j) - C_est_i(:,:,j) * R(:,:,i) * t;
        C_est_i(:,:,j) = C_est_i(:,:,j) * R(:,:,i);
    end
    
    if min(cost_i) < min(cost)
        index = find(cost_i == min(cost_i));
        C_est = C_est_i(:,:,index);
        t_est = t_est_i(:,index);
        cost = cost_i;
        flag = flag_i;
    end
    
    
end

if length(cost) < 1
    pause;
end

end

function r = rotx(t)
ct = cos(t);
st = sin(t);
r =    [1	0	0;
    0	ct	-st;
    0	st	ct];
end

function r = roty(t)
% roty: rotation about y-axi-
ct = cos(t);
st = sin(t);
r =    [ct	0	st;
    0	1	0;
    -st	0	ct];
end

function r = rotz(t)
% rotz: rotation about z-axis

ct = cos(t);
st = sin(t);
r = [ct	-st	0
    st	ct	0
    0	0	1];

end
