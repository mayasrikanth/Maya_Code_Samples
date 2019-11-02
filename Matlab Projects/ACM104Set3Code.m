%% (3) K-Means Clustering 
%%
struct = load("clustering_data.mat");
vectors = struct.x; 
clusters = zeros(1, 600); 
% Seeding 
rng(2016)
% Cluster representatives' indices 
r = randi([1, 600], 1, 3); 
 % Extracting representatives 
k_1 = vectors(r(1), 1:end);
k_2 = vectors(r(2), 1:end);
k_3 = vectors(r(3), 1:end);
% Storing representatives 
R = [k_1; k_2; k_3];
converge = false;
iter_counter = 0;
% Use for Plotting Original Data: 
vec_x = vectors(1:end,1);
vec_y = vectors(1:end,2);
% Contains values of p (objective function) 
all_p = [];
% Contains values of all cluster representatives
all_R = [];
while converge == false
    % Save plot if iter_counter = 1,5,or 10
    if iter_counter == 1 || iter_counter == 5 || iter_counter == 10
        % Original data 
        figure;
        %title("Cluster Image");
        scatter(vec_x, vec_y, 'black');
        hold on; 
        % Cluster 1
        %c_1 = clusters(clusters == 1);
        vec_x1 = vectors(clusters == 1, 1);
        vec_y1 = vectors(clusters == 1,2);
        scatter(vec_x1, vec_y1, 'red', 'filled');
        hold on;
        % Plotting Cluster 2
        scatter(vectors(clusters == 2, 1), vectors(clusters == 2, 2), 'green', 'filled');
        hold on; 
        % Plotting Cluster 3
        scatter(vectors(clusters == 3,1), vectors(clusters == 3, 2), 'blue', 'filled');
        %temp = strcat('Clusters', int2str(iter_counter)); 
        %title_format = 'Clusters%d';
        hold on;
        str_temp = sprintf('Clusters%d', iter_counter);
        %xlabel(str_temp);
        % Save the figure 
        %print(title, '-dpng');
    end 
    % Update c_prev
    c_prev = clusters; 
    p = 0;
    % for all vectors, update clustering  
    for k=1:600
       % UPDATE CLUSTERS
       currVec = vectors(k, 1:end);
       tempVec = zeros(1,3); 
       for i=1:3
           diff = currVec - R(i, 1:end);
           norm_diff = norm(diff);
           % recall these are 2-d
           res = (1/2) * norm_diff^2;
           tempVec(i) = res;
       end 
       % minimum cluster distance (objective function p)
       [m, idx] = min(tempVec);
       % Summing up the minimized objective function value
       % for each vector 
       p = p + m;  
       % storing cluster label for v_k
       clusters(k) = idx; 
    end 
    % Average p
    p_avg = p / 600;
    % Update vector of p values 
    all_p = [all_p; iter_counter + 1, p_avg];
    % UPDATE REPRESENTATIVES 
    idx_1 = find(clusters == 1); 
    idx_2 = find(clusters == 2);
    idx_3 = find(clusters == 3);
    % List of indices corresponding to vectors in cluster k 
    l_1 = length(idx_1);
    l_2 = length(idx_2);
    l_3 = length(idx_3);
    % For all vectors v_k s.t. c_k = 1, find average. 
    sum1 = [0 0];
    for i=1:l_1
        row_idx = idx_1(i);
        currvec = vectors(row_idx, 1:end);
        sum1 = sum1 + currvec;
    end
    % Store average R(1)
    avg1 = (1/l_1) * sum1; 
    R(1, 1:end) = avg1; % should be (x_avg, y_avg)
 
    sum2 = [0 0];
    for i=1:l_2
        row_idx = idx_2(i);
        currvec = vectors(row_idx, 1:end);
        sum2 = sum2 + currvec;
    end
    % Store average R(2)
    avg2 = (1/l_2) * sum2; 
    R(2, 1:end) = avg2; % should be (x_avg, y_avg)
    
    sum3 = [0 0];
    for i=1:l_3
        row_idx = idx_3(i);
        currvec = vectors(row_idx, 1:end);
        sum3 = sum3 + currvec;
    end
    % Store average R(3)
    avg3 = (1/l_3) * sum3; 
    R(3, 1:end) = avg3; % should be (x_avg, y_avg)
    
    % Updating Vector of cluster representatives 
    all_R = [all_R; R];
    % Testing for convergence: has the clustering changed? 
    if isequal(c_prev, clusters)
        converge = true; 
        break; 
    end 
    iter_counter = iter_counter + 1;
end 
% Plot the final clustering 
 % Original data 
figure;
scatter(vec_x, vec_y, 'black');
hold on; 
% Cluster 1
%c_1 = clusters(clusters == 1);
vec_x1 = vectors(clusters == 1, 1);
vec_y1 = vectors(clusters == 1,2);
scatter(vec_x1, vec_y1, 'red', 'filled');
hold on;
% Plotting Cluster 2
scatter(vectors(clusters == 2, 1), vectors(clusters == 2, 2), 'green', 'filled');
hold on; 
% Plotting Cluster 3
scatter(vectors(clusters == 3,1), vectors(clusters == 3, 2), 'blue', 'filled');
iter_counter = iter_counter + 1;
temp = strcat('Clusters', int2str(iter_counter)); 
title(temp);
% Save the figure 
%print(title, '-dpng');

% Plotting values of objective function p vs. iteration number 
figure;
scatter(all_p(:,1), all_p(:,2), 'green', 'filled');
temp2 = ('Objective Function Value vs. Iteration Number');
xlabel('Iteration Number');
ylabel('Objective Function Value');
title(temp2);
%title('Object Function vs. Iter Number');
%print('ObjectiveFunction', '-dpng');

% Minimum value achieved in last iteration? 
y_p = all_p(:,2);
min_p = y_p(length(y_p));
disp("Minimum value of objective function in last iteration: " + ...
num2str(min_p))

% Plotting trajectories of cluster representatives 
figure; 
scatter(vec_x, vec_y, 'black', 'filled');
hold on; 
scatter(all_p(:,1), all_p(:,2), 'blue', 'filled');
xlabel('X');
ylabel('Y');
%print('CenterTrajectories', '-dpng');

%% Comparing to Matlab's k-Means Implementation 
%%
matlab_idx = kmeans(vectors,3);
 % Original data 
figure;
scatter(vec_x, vec_y, 'black');
hold on; 
% Cluster 1
%c_1 = clusters(clusters == 1);
vec_x1 = vectors(matlab_idx == 1, 1);
vec_y1 = vectors(matlab_idx == 1, 2);
scatter(vec_x1, vec_y1, 'red', 'filled');
hold on;
% Plotting Cluster 2
scatter(vectors(matlab_idx == 2, 1), vectors(matlab_idx == 2, 2), 'green', 'filled');
hold on; 
% Plotting Cluster 3
scatter(vectors(matlab_idx == 3,1), vectors(matlab_idx == 3, 2), 'blue', 'filled');
title = "MatlabKMeans";
xlabel("Matlab X");
ylabel("Matlab Y");
% Save the figure 
print(title, '-dpng');
