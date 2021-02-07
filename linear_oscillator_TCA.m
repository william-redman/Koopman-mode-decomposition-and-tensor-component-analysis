%% Linear Oscilattor TCA
%-------------------------------------------------------------------------%
%   This script executes the linear oscillator examples from Mezic's
%   textbook (also similar to Example 13 in Budisic, Mohr, Mezic 2012) and 
%   compares the true Koopman mode decomposition (KMD) with that obtained
%   by using tensor component analysis (TCA). 
%
%   It uses TensorLab 3.0 (https://www.tensorlab.net/). 
%
%   Example #1 in "On Koopman Mode Decomposition and Tensor Component 
%   Analysis" W. T. Redman (2021). 
%   
%   Written by WTR 09/18/2020 // Last updated by WTR 01/07/2021
%-------------------------------------------------------------------------%
close all
%% Globals
c = 1;                                      % Damping term (note, it should be positive for damping)
A = [-c 1; -1 0];                          
dt = 0.001;                                 % Simulation time step
T_sim = 10000;                              % Number of time steps used in simulation
T = 1000;                                   % Number of time steps used for TCA and KMD (I found setting this not equal to T_sim led to slightly cleaner results)
n_cond = 20;                                % Number of initial conditions 
init_cond = cat(1, linspace(-1, 1, n_cond), zeros(1, n_cond));

%% Simulating the linear oscillator
x = zeros(2, T_sim, n_cond); 
x(:, 1, :) = init_cond; 

for tt = 2:T_sim
    x(:, tt, :) = A * squeeze(x(:, tt - 1, :)) * dt + squeeze(x(:, tt - 1, :)); 
end
x = x(:, (T_sim / T):(T_sim / T):T_sim, :); 

%% Finding the true KMD
%   This gives us three sets of vectors: lambda, phi, and v 
[v_nonabs, u] = eig(A);
[w, ~] = eig(A'); 
sign_v = sign(imag(v_nonabs)); 
sign_v(sign_v == 0) = 1; 
v = abs(v_nonabs);

phi = init_cond' * w; 
lambda = zeros(2, T_sim); 
lambda_nonreal = zeros(2, T_sim); 
lambda_nonreal(:, 1) = diag(exp(u * 0)); 
lambda(:, 1) = real(lambda_nonreal(:, 1));

for ii = 2:T_sim 
    lambda_nonreal(:, ii) = diag(exp(u * dt * (ii - 1))); 
    lambda(:, ii) = real(lambda_nonreal(:, ii)); 
end

lambda = lambda(:, (T_sim / T):(T_sim / T):T_sim);  

%% TCA modes 
%   This uses TensorLab 3.0 and follows the manual available on their site.  
R = 2;                                      % Number of TCA modes
model = struct; 
model.variables.a = randn(size(x, 1), R);
model.variables.b = randn(size(x, 2), R); 
model.variables.c = randn(size(x, 3), R); 
model.factors.A = {'a'}; 
model.factors.B = {'b'};
model.factors.C = {'c'}; 
model.factorizations.tensor.data = x;
model.factorizations.tensor.cpd = {'A', 'B', 'C'}; 
sdf_check(model, 'print'); 
sol = sdf_nls(model); 
[Uhat, output] = cpd(x, R); 

system_factors = sol.factors.A; 
time_factors = sol.factors.B; 
init_cond_factors = sol.factors.C; 

%% Finding optimal ordering
% The modes that TCA gives can be in a different order than the true KMD
% modes, so this find the corresponding modes (by maximum correlation).
% Made it slightly more complicated so it can handle when R > 2.

ordering = zeros(1, 2);
options = 1:R;

for ii = 1:2
    corr_KMD_mode_ii = zeros(1, length(options)); 
    for jj = 1:length(options)
        corr_KMD_mode_ii(jj) = corr(lambda(ii, :)', time_factors(:, options(jj))); 
    end
    [~, max_id] = max(corr_KMD_mode_ii); 
    ordering(ii) = options(max_id); 
    options(max_id) = []; 
    
end

%% Plotting results
for ii = 1:2
    % True modes
    figure %saving for bar plot later
    hold on
    
    figure
    plot(lambda(ii, :), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on
    
    figure
    plot(real(phi(:, ii)), 'ko', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on 
    
    % TCA modes
    [l_max, l_id] = max(abs(lambda(ii, :))); 
    [t_max, t_id] = max(abs(time_factors(:, ordering(ii)))); 
    time_scaling = sign(lambda(ii, l_id)) * sign(time_factors(t_id, ordering(ii))) * l_max / t_max; % Amount to rescale the time dependence
    
    [p_max, p_id] = max(abs(real(phi(:, ii)))); 
    [i_max, i_id] = max(abs(init_cond_factors(:, ordering(ii)))); 
    init_scaling = sign(real(phi(p_id, ii))) * sign(init_cond_factors(i_id, ordering(ii))) * p_max / i_max; % Amount to scale initial condition factors 
        
    figure(3 * (ii - 1) + 1)
    bar([1, 2], [v(:, ii), abs(system_factors(:, ordering(ii)) / (init_scaling * time_scaling))]); hold on 
    title(strcat({'Koopman modes / participation factors for mode'}, {' '}, {num2str(ii)}))
    legend('KMD', 'TCA'); 
    
    figure(3 * (ii - 1) + 2)
    plot(time_factors(:, ordering(ii)) * time_scaling, 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on
    title(strcat({'Eigenvalues / time dependence for mode'}, {' '}, {num2str(ii)}))
    legend('KMD', 'TCA'); 
    
    figure(3 * (ii - 1) + 3)
    plot(init_cond_factors(:, ordering(ii)) * init_scaling, 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on
    title(strcat({'Eigenfunctions / initial condition factors for mode'}, {' '}, {num2str(ii)}))    
    legend('KMD', 'TCA'); 
end

