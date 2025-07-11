function [S_E, I_E, S_I, I_I, eta_E, eta_I] = EI_dMFM(SC, dt, T, w, G, sigma, H_E, H_I, tau, varargin)

n = length(SC);

p = inputParser;
p.addParameter('J',0.15); % nA
p.addParameter('I_b',0.382);    % nA
p.addParameter('gamma',0.641);
p.addParameter('S_E_0',linspace(0,1,n));
parse(p,varargin{:}); 

J = p.Results.J;
I_b = p.Results.I_b;
gamma = p.Results.gamma; 


% Synaptic timescale
tau_e = tau(:,1);
tau_i = tau(:,2);

% Synaptic Weight
w_E = w(:,1);
w_I = w(:,2);
w_IE = w(:,3);
w_EI = w(:,4:3+29);
w_EE = w(:,33:end);

%% Simulation
tpre = ceil(10/dt);
tpost = ceil(T/dt);

% S_E = 0.164757 * ones([n, tpost+tpre]);
S_E = zeros([n, tpost+tpre]);
S_I = zeros([n, tpost+tpre]);
eta_E = S_E;
eta_I = S_I;
I_E = zeros([n, tpost+tpre]);
I_I = zeros([n, tpost+tpre]);

% EI_dMFM
for t=1:tpost+tpre-1
    eta_E(:,t) = sigma.*randn([n 1]).*sqrt(dt);
    eta_I(:,t) = sigma.*randn([n 1]).*sqrt(dt);
    I_E(:,t+1) = w_E.*I_b + w_EE'.*S_E(:,t) + G.*J.*SC*S_E(:,t) - w_IE.*S_I(:,t);
    I_I(:,t+1) = w_I.*I_b + w_EI'.*S_E(:,t) - S_I(:,t);
    % I_I(:,t+1) = w_i.*I_b + G.*J.*SC*S_E(:,t) - S_I(:,t);

    S_E(:,t+1) = S_E(:,t) + dt.*(-S_E(:,t)./tau_e + ...
     gamma.*(1-S_E(:,t)).*H_E(I_E(:,t+1))) + eta_E(:,t);
    S_I(:,t+1) = S_I(:,t) + dt.*(-S_I(:,t)./tau_i + ...
     H_I(I_I(:,t+1))) + + eta_I(:,t);
    
    % S_E(S_E(:,t+1)<0,t+1) = 0;
    % S_E(S_E(:,t+1)>1,t+1) = 1;
    % S_I(S_I(:,t+1)<0,t+1) = 0;
    % S_I(S_I(:,t+1)>1,t+1) = 1;
end

S_E = S_E(:, tpre+1:end);
I_E = I_E(:, tpre+1:end);
S_I = S_I(:, tpre+1:end);
I_I = I_I(:, tpre+1:end);

end

