%%
% Feb. 25 2021, He Zhang, fuyinzh@gmail.com 
% 
%   retrieve data from tracked feature data 
%

function [XXw, xxn] = get_data(filename)
    
    A = load(filename); 
    d = A(:,1); 
    
    XXw = [A(:,2), A(:,3), ones(size(d,1),1)]; 
    XXw = XXw.*d; 
    xxn = A(:,4:end); 

end
