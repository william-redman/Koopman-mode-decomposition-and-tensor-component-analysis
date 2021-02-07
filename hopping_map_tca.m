%% Hopping map TCA
%-------------------------------------------------------------------------%
%   This script executes the hopping map described in Ex. 1 of Mezic 2005
%   Nonlinear Dynamics and compares the true Koopman mode decomposition 
%   (KMD) with that obtained by using tensor component analysis (TCA).
%
%   It uses TensorLab 3.0 (https://www.tensorlab.net/). 
%
%   Example #2 in "On Koopman Mode Decomposition and Tensor Component 
%   Analysis" W. T. Redman (2021). 
%   
%   Written by WTR 09/07/2020 // Last updated by WTR 01/07/2021
%-------------------------------------------------------------------------%
%% Globals 
n_cond = 10;                                % Number of potential initial conditions
init_cond = sqrt(2) * ((-1):(1/n_cond):1);  
init_cond(abs(init_cond) > 1) = []; 
init_cond(init_cond == 0) = [];             % Generating initial conditions that are smaller than 1 (given the domain of the map) 
T = 30;                                     % Number of time steps used in simulation

%% Simulating the Hopping map 
traj = zeros(T + 1, length(init_cond)); 
traj(1, :) = init_cond; 

for tt = 2:(T + 1)
    for ii = 1:length(init_cond)
        x = traj(tt - 1, ii); 
        if x > 0.5
            traj(tt, ii) = 2 * (x - 0.5) - 1; 
        elseif x > 0
            traj(tt, ii) = 2 * x - 1;        
        elseif x < -0.5
            traj(tt, ii) = -2 * (x + 0.5);
        else
            traj(tt, ii) = -2 * x; 
        end
    end
end

X = reshape(traj, [1, T + 1, length(init_cond)]); 

%% TCA modes
%   This uses TensorLab 3.0 and follows the manual available on their site.  
R = 1;                                      % Number of TCA modes
model = struct; 
model.variables.a = randn(size(X, 1), R);
model.variables.b = randn(size(X, 2), R); 
model.variables.c = randn(size(X, 3), R); 
model.factors.A = {'a'}; 
model.factors.B = {'b'};
model.factors.C = {'c'}; 
model.factorizations.tensor.data = X;
model.factorizations.tensor.cpd = {'A', 'B', 'C'}; 
[Uhat, output] = cpd(X, R); 
sdf_check(model, 'print'); 
sol = sdf_nls(model); 

system_factors = sol.factors.A; 
time_factors = sol.factors.B; 
init_cond_factors = sol.factors.C; 
    
%% Plotting TCA vs KMD results
figure 
true_lambda = ones(1, T); 
true_lambda(2:2:end) = -1; 
true_lambda = true_lambda * sign(time_factors(1, 1)); 

plot(true_lambda, 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on 
plot(time_factors(:, 1) / mean(abs(time_factors(:, 1))), 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); 

title(strcat('Eigenvalues / time dependence for mode 1'))
legend('KMD', 'TCA'); 

figure
init_cond_mean = mean(abs(init_cond_factors(:, 1)));

plot(init_cond_mean .* sign(init_cond_factors(:, 1)), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on 
plot(init_cond_factors(:, 1), 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');  
title('Eigenfunctions / initial condition factors for mode 2')    
legend('KMD', 'TCA'); 


