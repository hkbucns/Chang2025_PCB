% It is built for linear reconstruction based on cov and partial cov.
% The idea of this function is to regress y, tau and C based on the
% knowledge of Symmetric SC, Estimated Jacobian Matrix and Stable Points.

% Firstly, we will build Matrix Mat_cont for regression, that is, to comput
% Matrix_cont * y = S_vec, where Matrix_cont, y and S_vec stand for the
% Symmetric Constraints, Fixed-points related Heterogeneity and Vectorized
% Symmetric SC.
% After getting y, we can reveal timescale Tau and 
% Asymmetric Structural Connectivity C.

% Input:
% SC: prior knowledge from Structural Connectivity (Symmetric);
% W: Reconstructed Jacobian Matrix(Have same size as S).

% Output:
% H: Reconstructed hierarchical heterogeneity,
% C_recon: Reconstructed connectivity.

% Inter Parameters:
% S_vec: reshaped S, easier for estimation;
% W_trans: Transpose of W;
% Mat_cont: Temporary Matrix for defining R;
% R: Concatenated Matrix for revealing H;
% S_index: non-zero element indeices;
% y_st: 1/(gamma*G*J*(1-S_star)*h_i), where h_i is the slope of H_i at x_i;

function [y_st, C_recon] = RevealHHetero2(SC, W)
N = size(W,1);
SC_vec = 2*reshape(SC,N^2,1);
W_trans = W';
Mat_cont = zeros(N);
Mat_cont(:,1) = W_trans(:,1);
Mat_cont = Mat_cont+diag(W(:,1));
R = Mat_cont;
for k = 2:N
    Mat_cont = zeros(N);
    Mat_cont(:,k) = W_trans(:,k);
    Mat_cont = Mat_cont+diag(W(:,k));
    R = cat(1,R,Mat_cont);
end
% Eliminating Diagonal Elements
SC_index = (SC_vec~=0);
SC_vec = SC_vec(SC_index);
R1 = R(SC_index,:); 

% Use lsqr for iterative solution
maxit = 500; % You can adjust this value
tol = 1e-6; % You can adjust this value
y_st = lsqr(R1, SC_vec, tol, maxit);
C_recon = W .* (y_st * ones(1, N));
end
