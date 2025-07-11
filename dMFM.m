function [S,eta] = dMFM(SC, dt, T, w, I, G, sigma, varargin)
%   SC is the Structural Connectivity
%   dt is the time step, T is total running seconds
%   w is Self-recurrent Strength，I is External Input，
%   G is Global Coupling, sigma is Noise Strength

p = inputParser;            
p.addParameter('J',0.2609); % nA
p.addParameter('tau_s',0.1);    % s 
p.addParameter('gamma_s',0.641);
parse(p,varargin{:}); 

J = p.Results.J;
tau_s = p.Results.tau_s;
% tau_s = linspace(0.05,0.15,29)';
gamma_s = p.Results.gamma_s; 

a = 270; 
b = 108; % Hz
d = 0.154; % s
H = @(x)dMFM_H(x,a,b,d);

% pre run 10s.
n = length(SC); 
tpre = ceil(10/dt);  
if dt <= 0.001
    tpre = ceil(1/dt);  
end
tpost = ceil(T/dt); 

S = ones([n, tpost+tpre])/10;
eta = S;

for t=1:length(S)-1
     eta(:,t+1) = sigma.*randn([n 1]).*sqrt(dt);
     x = w.*J.*S(:,t) + G.*J.*SC*S(:,t) + I;
     S(:,t+1) = S(:,t) + dt.*(-S(:,t)./tau_s + ...
         gamma_s.*(1-S(:,t)).*H(x)) + eta(:,t+1);
     S(S(:,t+1)<0,t+1) = 0;
     S(S(:,t+1)>1,t+1) = 1;
end

S = S(:, tpre+1:end);
eta = eta(:,tpre+1:end);

end

