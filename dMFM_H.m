function [H] = dMFM_H(x, a, b, d)
%   Activation Function
%   https://www.jneurosci.org/content/33/27/11239
H = (a*x - b) ./ (1 - exp(-d*(a*x-b)));
end

