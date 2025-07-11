function [S,eta] = dMFM_eq(SC, dt, T, tau, b, G, sigma, varargin)


p = inputParser;           
p.addParameter('J',0.2609); % nA
p.addParameter('gamma_s',0.641);
parse(p,varargin{:}); 

J = p.Results.J;
gamma_s = p.Results.gamma_s; 

w = 0.5327;
I = 0.3266;

a = 270; %n/C (/nC?)
d = 0.154; % s
H = @(x)dMFM_H(x,a,b,d);

n = length(SC); 
tpre = ceil(10/dt);  
if dt <= 0.001
    tpre = ceil(1/dt); 
end
tpost = ceil(T/dt); 

S = ones([n, tpost+tpre])/10;
eta = zeros([n, tpost+tpre]);

for t=1:length(S)-1
     eta(:,t+1) = eta(:,t+1) + sigma.*randn([n 1]).*sqrt(dt);
     x_val = w.*J.*S(:,t) + G.*J.*SC*S(:,t) + I;
     H_vals = H(x_val);
     S(:,t+1) = S(:,t) + dt.*(-S(:,t)./tau + ...
         gamma_s.*(1-S(:,t)).*H_vals) + eta(:,t+1);
     S(S(:,t+1)<0,t+1) = 0;
     S(S(:,t+1)>1,t+1) = 1;
end

S = S(:, tpre+1:end);

end

