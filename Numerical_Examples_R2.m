%% Numerical examples: R = 2
%-------------------------------------------------------------------------%
%   This script generates the data used in Fig. 1 of "On Koopman Mode 
%   Decomposition and Tensor Component Analysis" Redman (2021). 
%
%   The Exact DMD modes were computed using Algorithm 2 and Appendix 1 of 
%   Tu et al. 2014. The TCA modes were computed using TensorLab 3.0. 
%
%   For more details see Sec. V of Redman 2021. Feel free to send any 
%   questions or comments to wredman@ucsb.edu. 
%
%   Written by WTR 03/06/2021 // Last updated by WTR 03/11/2021
%-------------------------------------------------------------------------%
%% Generating the data
dim = [2, 101, 10];                                                        % size of the data tensors 
n_modes = 2;                                                               % number of sources underlying the data

v_bar = rand(dim(1), n_modes);      
v_bar(:, 1) = v_bar(:, 1) / norm(v_bar(:, 1))^2;                           % making the modes orthonormal
v_bar(:, 2) = v_bar(:, 2) - dot(v_bar(:, 1), v_bar(:, 2)) * v_bar(:, 1) / norm(v_bar(:, 1))^2; 
v_bar(:, 2) = v_bar(:, 2) / norm(v_bar(:, 2))^2;

alpha = 2 * rand(1, n_modes) - 1; 
beta = 2 * rand(1, n_modes) - 1; 

true_mags = exp(2 * pi * 1i * alpha/10 + beta/10);                         % sorting the true modes by size of their imaginary component
[~, ids] = sort(imag(true_mags)); 

alpha = alpha(ids); 
beta = beta(ids); 

S_bar = zeros(n_modes, dim(2));                                            % generating the time evolution of each mode
for ii = 1:n_modes
    S_bar(ii, :) = (exp(2 * pi * 1i * alpha(ii)/10 + beta(ii)/10)).^(0:(dim(2) - 1)); 
end

gamma = 3 * rand(1, n_modes);                                              % generating the dependence on initial condition
Phi_bar = zeros(dim(3), n_modes); 
x = linspace(0.1, 1, dim(3)); 
for ii = 1:n_modes
    Phi_bar(:, ii) = x.^gamma(ii);
end

D = zeros(dim);                                                            % constructing the full data tensorr
for ii = 1:dim(3)
    for jj = 1:n_modes
        D(:, :, ii) = D(:, :, ii) + v_bar(:, jj) * S_bar(jj, :) * Phi_bar(ii, jj);
    end
end

%% Exact DMD
X = D(:, 1:(end - 1), :); 
Y = D(:, 2:end, :); 

v_DMD = zeros(dim(1), n_modes, dim(3)); 
S_DMD = zeros(dim(1), dim(2) - 1, dim(3)); 
Phi_DMD = zeros(n_modes, dim(3)); 

for ii = 1:dim(3) 
    [U, S_svd, V] = svd(X(:, :, ii), 'econ'); 
    A_tilde = U' * Y(:, :, ii) * V * inv(S_svd); 
    
    [w,u] = eig(A_tilde);
     
    if ((imag(u(1,1)) > imag(u(2,2))) + (alpha(1) > alpha(2))) == 1        % sorting the DMD modes to match the true modes
        w = w(:, [2, 1]); 
        u = [u(2, 2), 0; 0, u(1,1)]; 
    end
    
    S_DMD(:, :, ii) = diag(u).^(1:(dim(2) - 1)); 
    
    for jj = 1:n_modes
        if abs(u(jj, jj)) > 10^(-10)
            Koop_mode = (1 ./ u(jj, jj)) * Y(:, :, ii) * V * inv(S_svd) * w(:, jj); 
            v_DMD(:, jj, ii) = Koop_mode * v_bar(1, jj) / Koop_mode(1);            
        end
    end
    
    Phi_DMD(:, ii) = inv(u) * pinv(v_DMD(:, :, ii)) * Y(:, 1, ii); 
    
end

%% TCA
R = n_modes;                                                               % number of TCA modes
model = struct; 
model.variables.a = randn(size(Y, 1), R) + 1i * randn(size(Y, 1), R); 
model.variables.b = randn(size(Y, 2), R) + 1i * randn(size(Y, 2), R);  
model.variables.c = randn(size(Y, 3), R) + 1i * randn(size(Y, 3), R); 
model.factors.A = {'a'}; 
model.factors.B = {'b'};
model.factors.C = {'c'}; 
model.factorizations.tensor.data = Y;
model.factorizations.tensor.cpd = {'A', 'B', 'C'}; 
sdf_check(model, 'print'); 
sol = sdf_nls(model); 
[Uhat, output] = cpd(Y, R); 
relerr = output.Algorithm.relerr                                           % relative error in TCA decomposition (as compared to the true Y)

A = sol.factors.A; 
B = sol.factors.B; 
C = sol.factors.C;  

% Ordering the modes 
CC = corr(real(Phi_bar), real(C));
[~, max_id] = max(abs(CC(:))); 
if max_id ~= 1 && max_id ~= 4
    ordering = [2,1];
    A = A(:, ordering); 
    B = B(:, ordering); 
    C = C(:, ordering);     
end

% Scaling the modes       
for ii = 1:R
    B_scaling = S_bar(ii, 2) / B(1, ii); 
    C_scaling = Phi_bar(1, ii) / C(1, ii); 
    B(:, ii) = B(:, ii) * B_scaling;
    C(:, ii) = C_scaling * C(:, ii); 
    A(:, ii) = A(:, ii) / (B_scaling * C_scaling); 
end

%% Plotting
% The TCA and DMD modes are offset from each other. To see them ontop of
% each other, uncomment lines 148 and 149. 
for ii = 1:n_modes
    figure
    bar([1, 2], [v_bar(:, ii), real(mean(v_DMD(:, ii, :), 3)), real(A(:, ii))]); 
    title(strcat({'Mode ',}, {num2str(ii)}, ' V'));
    
    figure
    plot(real(S_bar(ii, 2:end)), 'k-', 'LineWidth', 2); hold on
    plot(2:2:size(S_DMD, 2), real(mean(S_DMD(ii, 2:2:end, :), 3)), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on 
    plot(3:2:size(B, 1), real(B(3:2:end, ii)), 'o','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); 
    title(strcat({'Mode ',}, {num2str(ii)}, ' Lambda'));
    
    figure
    shift = 0.01;
    plot([(0.1 - 2 * shift):0.01:(1 + 2 * shift)], [(0.1 - 2 * shift):0.01:(1 + 2 * shift)].^gamma(ii) , 'k-'); hold on 
    true_modes = x.^gamma(ii);   
    TCA_error = reshape(C(:, ii), size(true_modes)) - true_modes;
    DMD_error = Phi_DMD(ii, :) - true_modes; 
    TCA_x = x + shift; 
    DMD_x = x - shift; 
    TCA_y = TCA_error .* (TCA_x ./ x).^gamma(ii) + TCA_x.^gamma(ii); 
    DMD_y = DMD_x.^gamma(ii) - DMD_error .* (DMD_x ./ x).^gamma(ii); 
    plot(TCA_x, TCA_y, 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    plot(DMD_x, DMD_y, 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on
%     plot(0.1:0.1:1, abs(C(:, ii)), 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%     plot(0.1:0.1:1, abs(Phi_DMD(ii, :)),'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); hold on
    title(strcat({'Mode ',}, {num2str(ii)}, ' Phi'));
    axis([0.1 - 2 * shift, (1 + 2 * shift), 0, 1.05])
end



