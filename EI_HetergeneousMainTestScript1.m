%% Parameters setting

load('SC.mat');
SC = fln;
SC = SC./max(SC,[],'all');
N = length(SC);
T = 5000; %s, ;
dt = 0.001; %1ms, ;
t = linspace(0,T,T/dt);

tau = [0.1 , 0.01]; % E and I
gamma = 0.641;
G = 1.1; 
sigma = 0.01; % nA, the noise amplitud.
w_E = 1; % Scale I_0
w_I = 0.7; % Scale I_0
J = 0.15; %J_NMDA
I_0 = 0.382;
w_IE = 1; %nA, J_i.         % Can be adjusted to set firing rate ~ 3Hz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Heterogeneity
w_EE = 0.21*linspace(0.6,1,N);
% Local excitatory recurrent, w_plus*J_NMDA = 0.21;

w_EI = 0.15*linspace(0.6,1,N);
% Used to be settled as heterogeneous.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = [w_E, w_I, w_IE, w_EI, w_EE];
%% Setting f-I function H and dH

a_E = 310; % nC
b_E = 125; % Hz
d_E = 0.16; % s

a_I = 615;
b_I = 177;
d_I = 0.087;

H_E = @(x)dMFM_H(x,a_E,b_E,d_E); % f-I Curve for Excitatory Population
H_I = @(x)dMFM_H(x,a_I,b_I,d_I); % f-I Curve for Inhibitory Population

dH_E = @(x) - 310./(exp(20 - (248.*x)./5) - 1) - ...
    (248.*exp(20 - (248.*x)./5).*(310.*x - 125))./(5*(exp(20 - (248.*x)./5) - 1).^2);
dH_I = @(x) - 615./(exp(15399/1000 - (10701.*x)./200) - 1) - ...
    (10701.*exp(15399/1000 - (10701.*x)./200).*(615.*x - 177))./(200*(exp(15399/1000 - (10701.*x)./200) - 1).^2);

%% Test the matching of H and dH

% t = linspace(0,1,1000);
% H_plot = H(t);
% dH_num = (H_plot(2:end)-H_plot(1:end-1))/0.001;
% dH_ana = dH(t(2:end));
% figure(1)
% plot(t(1:end-1),dH_num,t(1:end-1),dH_ana);
% title('numerical dH and analytical dH');
% figure(2)
% plot(t,H_plot);
% title('numerical H')
% % Good derivative

%% Reduced Wong-Wang Model Simulation

[S_E, I_E, S_I, I_I] = EI_dMFM(SC, dt, T, w, G, sigma, H_E, H_I, tau);
S_E_star = mean(S_E,2);
S_I_star = mean(S_I,2);


x_E_star = w_EE'.*S_E_star + G.*J.*SC*S_E_star - w_IE.*S_I_star + w_E*I_0;
x_I_star = w_EI'.*S_E_star - S_I_star + w_I*I_0;
    
%% Jacobian Matrix Calculation

Jacob_EE = zeros(N);
Jacob_EI = zeros(N);
Jacob_IE = zeros(N);
Jacob_II = zeros(N);
for i = 1:N
    for j = 1:N
        if i == j
            Jacob_EE(i,j) = -1/tau(1)-gamma*H_E(x_E_star(i))+w_EE(i)*gamma*(1-S_E_star(i))*dH_E(x_E_star(i));
            Jacob_EI(i,j) = w_EI(i)*dH_I(x_I_star(i));
            Jacob_IE(i,j) = -w_IE*(1-S_E_star(i))*gamma*dH_E(x_E_star(i));
            Jacob_II(i,j) = -1/tau(2) - dH_I(x_I_star(i));
        else
            Jacob_EE(i,j) = gamma*G*J*(1-S_E_star(i))*SC(i,j)*dH_E(x_E_star(i));
        end
    end
end
Jacob = [Jacob_EE Jacob_IE
    Jacob_EI Jacob_II];

%% Linear Simulation & Validation
% 
% tpre = ceil(10/dt);
% tpost = ceil(T/dt);
% % S_lin = [0.164757 * ones([N, tpost+tpre]);zeros([N, tpost+tpre])];
% S_lin = [zeros([N, tpost+tpre]);zeros([N, tpost+tpre])];
% S_star = [S_E_star;S_I_star];
% eta = [eta_E;eta_I];
% for t=1:length(S_lin)-1
%      S_lin(:,t+1) = S_lin(:,t) + dt.*(Jacob*(S_lin(:,t)-S_star)) + eta(:,t);
%      S_lin(S_lin(:,t+1)<0,t+1) = 0;
%      S_lin(S_lin(:,t+1)>1,t+1) = 1;
% end
% S_lin_plot = S_lin(:, tpre+1:end);
% clear eta tpre tpost t 
% 
% % Showing Simulations
% figure(3)
% t = linspace(0,T,T/dt);
% plot(t,S_lin_plot)
% set(gca,'box','off');
% title('Linear Simulation')

% figure(4)
% plot(t,S_E)
% set(gca,'box','off');
% title('Nonlinear Simulation E')
% figure(3)
% plot(t,S_I)
% set(gca,'box','off');
% title('Nonlinear Simulation I')

%% Linear Reconstruct Jacobian Matrix Jacob_eff with S_E
% S = [S_E;S_I];
% Jacob_est = LinearReconst(S,dt); % Estimating Jacobian Matrix
Jacob_est = LinearReconst(S_E,dt);
% Jacob_est_EE = Jacob_est(1:N,1:N);

%% Reconstruct spatial properties with Jacobian Matrix
SC_sym = (SC + SC')/2;
% SC_sym = SC_sym.*(SC_sym>1e-3);
% [y_st,w_recon,C_recon,R1] = RevealHHetero1(SC_sym,Jacob_est_EE,S_E_star,tau(1),gamma,G, J);
[y_st, C_recon] = RevealHHetero2(SC_sym,Jacob_est); % A compact version of spatial reconstruction
C_recon = C_recon-diag(diag(C_recon));

%% Evaluation of Jacobian Matrix J and SC Estimation
[SSE_J,Corr_J,Corr_nonzero_J] = EstimationJacobianPlotting(Jacob_EE,Jacob_est);
[SSE_SC,Corr_SC,Corr_nonzero_SC] = EstimationMatrixPlotting(SC,C_recon);

%% Effective Heterogeneity
figure; scatter(gamma*G*J*(1-S_E_star).*dH_E(x_E_star),1./y_st);