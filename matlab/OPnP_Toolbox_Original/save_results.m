%%
% Feb. 17, 2021, write the results to files 

function save_results(xs, method_list)
    
    fmean = fopen('test_results/mean.txt', 'w'); 
    % fmed = fopen('test_results/median.txt', 'w');
    fprintf(fmean, 'number of points [mean_t, std_t, mean_r, std_r] \n');
    for i= 1:length(method_list)
        fprintf(fmean, ' %s ', method_list(i).name); 
    end
    fprintf(fmean, '\n');
    
    for j =1:length(xs)
        fprintf(fmean, '%d ', xs(j));
        for i= 1:length(method_list)
            fprintf(fmean, ' %4.8f %4.8f %4.8f %4.8f ', method_list(i).mean_t(j), ...
                method_list(i).std_t(j), method_list(i).mean_r(j), ... 
                method_list(i).std_r(j)); 
        end
        fprintf(fmean, '\n'); 
    end
    fclose(fmean);
end