% load data
T_covid = readtable('matlab/data_covid.csv');
T_non_covid = readtable('matlab/data_non_covid.csv');
% estimate parameters (constant)
d = T_covid(:, [2 4 5 6]);
d_n = T_non_covid(:, [2 4 5 6]);
guess = [60, 0.5];

% set options
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',3000);
options.MaxIterations = 3000;
options.FunctionTolerance = 1e-10;
options.StepTolerance = 1e-10;
lb = [-Inf, 1e-6];

%estimate
estimates_non_covid = lsqnonlin(@SSR, guess, lb, [], options, d_n{:,:})
estimates_covid = lsqnonlin(@SSR, guess, lb, [], options, d{:,:})

% create function for SSR for mean-reversion model
function f = SSR(c, x)
    dif = x(:,2) - x(:,1);
    Fhat = c(:,1) + dif.*exp(-1*c(:,2).*x(:,4));
    f = (x(:,1) - Fhat).^2;
end