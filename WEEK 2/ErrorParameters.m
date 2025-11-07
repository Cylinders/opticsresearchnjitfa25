N = 1:1:150;
RMSE = zeros(size(N));
R2   = zeros(size(N));
MAE  = zeros(size(N));

for k = 1:length(N)
    [RMSE(k), R2(k), MAE(k)] = dopplererror(N(k));
end




% Create a table
T = table(N', RMSE', R2', MAE', ...
    'VariableNames', {'N', 'RMSE', 'R2', 'MAE'});

% Display the table
disp(T);

% (Optional) Save table to a file, e.g., CSV
writetable(T, 'doppler_error_metrics.csv');
