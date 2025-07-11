% It is built for linear reconstruction based on cov and partial cov.

% Input:
% x: time series
% dt: time step

% Output:
% L: series of reconstructed matrix, size N,N

function [L] = LinearReconst(x1,dt)
rate = 1;
if rate == 1
    x = x1;
else
    x = downsample(x1',rate)';
    dt = dt*rate;
end
[N,M]=size(x);
L = zeros(N);
dx = (x(:,2:M)-x(:,1:M-1))./dt;
C2 = cov([dx; x(:,1:M-1)]'); % Calculating Partial Cov
C1 = C2(1:N,N+1:2*N);
C = cov(x(:,1:M-1)');
L(:,:) = C1/C;
end