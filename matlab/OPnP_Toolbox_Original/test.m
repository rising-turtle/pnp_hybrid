

%% Camera Ci's pose 
Ri = eye(3); 
ti = [0, 0, 0]'; 

%% Camera Cj's true pose 
Rj = eye(3); 
tj = [0, 0, 5]'; 

%% Camera Cj's predicted pose, there is 4 meters error along z-axis 
Rj_hat = eye(3); 
tj_hat = [0, 0, 1]'; 

%% feature k [7, 8, 9] in Ci's coordinate system 
Pti = [5, 3, 10]'; 

%% measurement in Cj's
% omit rotation since they are identity 
% Rji = Rj'*Ri; 
% tji = Rj'*(ti-tj); 
% Ptj = Rji*Pti + tji; 
Ptj = Pti - tj; 

%% measurement zj_bar has no noise 
zj_bar = Ptj(3);

%% project on Cj's normalized plane 
ptj = [Ptj(1)/zj_bar, Ptj(2)/zj_bar]'; 

%% prediction on Cj's 
Ptj_hat = Pti - tj_hat; 
zj_predict = Ptj_hat(3); 

%% if divide zj_bar, the residual is zero  
ptj_hat = [Ptj_hat(1)/zj_bar, Ptj_hat(2)/zj_bar]'; 
fprintf('If divide by zj_bar, the pose error along Z-axis cannot be computed since the residual is: ')
residual_1 = ptj - ptj_hat

%% if divide zj_predict, the residual is not zero, which can be used to correct the pose error 
ptj_hat = [Ptj_hat(1)/zj_predict, Ptj_hat(2)/zj_predict]'; 
fprintf('If divide by zj_predict, the pose error along Z-axis can be computed by minimizing the residual: ')
residual_2 = ptj - ptj_hat 



