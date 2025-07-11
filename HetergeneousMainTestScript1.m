%% Parameters setting

load('SC.mat');
SC = fln;
N = length(SC);
T = 50000; %s;
dt = 0.01; %10ms;
w = linspace(0.8, 1.3, N)' * 0.5; % Local excitatory recurrent;
I = linspace(1, 1.167, N)' * 0.3; % nA, the overall effective external input;

G = 0.65;
sigma = 0.01; % nA, the noise amplitud.
J = 0.2609;
tau = 0.1;
gamma = 0.641;

%% Setting f-I function H and dH
a = 270; %n/C (/nC?)
b = 108; % Hz
d = 0.154; % s
H = @(x)dMFM_H(x,a,b,d); % f-I Curve
dH = @(x) - 270./(exp(2079/125 - (2079*x)./50) - 1) -...
     (2079.*exp(2079/125 - (2079*x)./50).*(270*x - 108))./(50*(exp(2079/125 - (2079*x)./50) - 1).^2);

%% Reduced Wong-Wang Model (Model A) Simulation

[S,eta] = dMFM(SC, dt, T, w, I, G, sigma);
S_star = mean(S,2);
x_star = w.*J.*S_star+G.*J.*SC*S_star+I;

dH_val = dH(x_star); % reveal the value of dH/dt at x_star

%% Jacobian Matrix Calculation

Jacob = zeros(N);
y_ana_inv = zeros(N,1);
for i = 1:N
    y_ana_inv(i) = gamma*G*J*(1-S_star(i))*dH(x_star(i));
    for j = 1:N
        if i == j
            Jacob(i,j) = -1/(tau*(1-S_star(i)))+w(i)*gamma*J*(1-S_star(i))*dH(x_star(i));
        else
            Jacob(i,j) = gamma*G*J*(1-S_star(i))*SC(i,j)*dH(x_star(i));
        end
    end
end
clear i j 

%% Linear Reconstruct Jacobian Matrix with S
Jacob_est = LinearReconst(S,dt); % Estimating Jacobian Matrix

%% Reconstruct spatial properties with Jacobian Matrix
% y_st is the inverse of effective heterogeneity.

SC_sym = (SC + SC')/2;
[y_st,w_recon,C_recon] = RevealHHetero1(SC_sym,Jacob_est,S_star,tau,gamma,G, J);
C_recon = C_recon-diag(diag(C_recon));

%% Evaluation of Jacobian Matrix J and SC Estimation

[SSE_J,Corr_J,Corr_nonzero_J] = EstimationJacobianPlotting(Jacob,Jacob_est);
[SSE_SC,Corr_SC,Corr_nonzero_SC] = EstimationMatrixPlotting(SC,C_recon);

%% Evaluation of Partial Firing Rate dH/dx Estimation

dH_st = 1./(gamma*G*J.*(1-S_star).*y_st);

figure(5)
scatter(dH_val,dH_st);
set(gca,'box','off');
xlabel('Ground Truth dH/dx_i^*');
ylabel('Estimated dH/dx_i^*');
xticks([0 300]);
yticks([0 300]);
alpha(0.8);
saveas(gcf,'hPerformance.png');

%% Evaluation of Local excitatory recurrent w Estimation

figure(6)
scatter(w,w_recon);
set(gca,'box','off');
hold on
plot([0.4,max(w)],[0.4,max(w)],'--','Color','k');
hold off
xlabel('Ground Truth w_i');
ylabel('Estimated w_i');
xticks([0.4 0.65]);
yticks([0.4 0.65]);
alpha(0.8);
saveas(gcf,'wPerformance.png');

%% Evaluation of Input Current x* Estimation

x_st = zeros(size(dH_st));
for i = 1:length(dH_st)
    x_initial_guess = 0.5;
    x_solution = fsolve(@(x) dH(x) - dH_st(i), x_initial_guess);
    x_st(i) = x_solution;
end

% Plotting Partial f_I Curve Mapping relationship 
figure(7)
scatter(x_st,dH_st);
hold on
xLine = linspace(min(x_st), max(x_st), 100);
yLine = dH(xLine);
plot(xLine, yLine, 'r-');
hold off
set(gca,'box','off');
xlabel('Estimated x_i^*');
ylabel('Estimated dH/dx_i^*');
xticks([0 0.55]);
yticks([0 270]);
alpha(0.8);
saveas(gcf,'PFIMapping.png');

% Plotting Estimated Input Current x* performance
figure(8)
scatter(x_star,x_st);
set(gca,'box','off');
xlabel('Ground Truth x_i^*');
ylabel('Estimated x_i^*');
xticks([0.3 0.55]);
yticks([0.3 0.55]);
alpha(0.8);
saveas(gcf,'x_starPerformance.png');

% Plotting f-I Curve Mapping relationship 
figure(9)
scatter(x_st,H(x_st));
hold on
xLine = linspace(min(x_st), max(x_st), 100);
yLine = H(xLine);
plot(xLine, yLine, 'r-');
hold off
set(gca,'box','off');
xlabel('Estimated x_i^*');
ylabel('Firing rate H');
xticks([0 0.55]);
yticks([0 3 10]);
alpha(0.8);
saveas(gcf,'FIMapping.png');

%% Evaluation of Outer Input I Estimation

I_recon = x_st - w_recon.*J.*S_star - G.*J.*C_recon*S_star;

figure(10)
scatter(I,I_recon);
hold on
plot([0.3,max(I)],[0.3,max(I)],'--','Color','k');
hold off
set(gca,'box','off');
xlabel('Ground Truth I_i');
ylabel('Estimated I_i');
xticks([0.3 0.35]);
yticks([0.3 0.35]);
alpha(0.8);
saveas(gcf,'IPerformance.png');

%% Testing the error contribution of components from estimated w, I

% Plotting Relative Error of w, I with Stable points S*
figure(11)
scatter(S_star,abs(w-w_recon)./w);
set(gca,'box','off');
xlabel('Stable Gating Variable S^*');
ylabel('Relative Error of w');
xticks([0 0.8]);
yticks([0 0.1 0.2]);
saveas(gcf,'ErrorW.png');

figure(12)
scatter(S_star,abs(I-I_recon)./I);
set(gca,'box','off');
xlabel('Stable Gating Variable S^*');
ylabel('Relative Error of I');
xticks([0 0.8]);
yticks([0 0.05 0.15]);
saveas(gcf,'ErrorI.png');

close all


%% heatmap
figure(13)
h1 = heatmap(SC,'Colormap',hot,'GridVisible','off');

figure(14)
h2 = heatmap(C_recon,'Colormap',hot,'GridVisible','off');

%% Replacement or tau & b from w & I

tau_eq_inv = 1/tau-(w_recon-0.5327).*(1-S_star)./y_st./G;

dH_sym = @(b, x) (2079*exp((77*b)/500 -...
    (2079*x)/50)*(b - 270*x))/(50*(exp((77*b)/500 - (2079*x)/50) - 1)^2) - 270/(exp((77*b)/500 - (2079*x)/50) - 1);

x_star_eq = 0.5327*J*S_star+G*J*C_recon*S_star+0.3266;

b_fit = zeros(1,N);
options = optimoptions('lsqcurvefit', 'Display', 'iter', 'MaxIterations', 5000);
b0 = 108;
for i = 1:N
    b_fit(i) = lsqcurvefit(dH_sym, b0, x_star_eq(i), dH_st(i), [], [], options);
end

% figure(16)
% plot(b_fit,'o');xlabel('ROIs');ylabel('Reconstructed b_i');set(gca,'box','off');

%% Illustration of replaced tau & b across ROIs

b_fit_an = b-a*J.*(w_recon-0.5327).*S_star-a.*(I_recon-0.3266);

figure(15)
Index = 1:1:N;
t2 = tiledlayout(1,1);
ax1 = axes(t2);
scatter(ax1,Index,1./tau_eq_inv,'MarkerEdgeColor',[199,35,54]./256,'MarkerFaceColor',[199,35,54]./256);
% alpha(ax1,0.8);
% legend('\tau_i')

ax2 = axes(t2);
scatter(ax2,Index,b_fit,'MarkerEdgeColor',[33,49,80]./256,'MarkerFaceColor',[33,49,80]./256)
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
hold on
hold off
% alpha(ax2,0.8);

xlabel(ax1,'ROIs');
ylabel(ax1,'Replaced Features');

xticks(ax1,[]);
yticks(ax1,[]);
xticks(ax2,[]);
yticks(ax2,[]);

% figure(15)
% scatter(b_fit_an,b_fit,'MarkerEdgeColor',[33,49,80]./256,'MarkerFaceColor',[33,49,80]./256);
% xlabel('Analytical b_i');
% ylabel('Fitted b_i');
% set(gca,'box','off');
% hold on
% plot([min(b_fit),max(b_fit)],[min(b_fit),max(b_fit)],'--','Color','k');
% plot(b,b,'^','Color','b')
% hold off
% xticks([96 108 116]);
% yticks([96 108 116]);

%% Bar Chart of tau and b

% figure(16)
% yyaxis left
% bar(tau_eq,'FaceColor',[199,35,54]./256,'EdgeColor',[199,35,54]./256);
% ylim([0.09,0.12]);
% alpha(0.8);
% 
% yyaxis right
% bar(b_fit,'FaceColor',[33,49,80]./256,'EdgeColor',[33,49,80]./256);
% ylim([95 120]);
% alpha(0.8);
% set(gca,'box','off')

%% Examination of Replacement Performance

tau_eq = 1./tau_eq_inv;
[S2] = dMFM_eq(SC, dt, T, tau_eq, b_fit', G, sigma);
S_star_eq = mean(S2,2);
FC1 = corr(S');
FC2 = corr(S2');

%% Equilivant Jacobian Matrix from tau & b

Jacob_eq = zeros(N);
x_star_eq_new = 0.5327*J*S_star_eq+G*J*C_recon*S_star_eq+0.3266;
for i = 1:N
    for j = 1:N
        if i == j
            Jacob_eq(i,j) = -1/(tau_eq(i)*(1-S_star_eq(i)))+0.5327*gamma*J*(1-S_star_eq(i))*dH_sym(b_fit(i),x_star_eq(i));
        else
            Jacob_eq(i,j) = gamma*G*J*(1-S_star_eq(i))*SC(i,j)*dH_sym(b_fit(i),x_star_eq(i));
        end
    end
end
clear i j 


%% Plotting Replaced Dynamics

figure(1)
scatter(S_star,S_star_eq,'MarkerEdgeColor',[33,49,80]./256,'MarkerFaceColor',[33,49,80]./256);
hold on
plot([0,0.7],[0,0.7],'--','Color','k');
hold off
xlabel('Original Dynamics S^*');
ylabel('Replaced Dynamics S^*');
xticks([0 0.6]);
yticks([0 0.6]);

figure(2)
heatmap(FC1-diag(diag(FC1)),'Colormap',hot,'GridVisible','off','ColorLimits',[0 0.7]);
figure(3)
heatmap(FC2-diag(diag(FC2)),'Colormap',hot,'GridVisible','off','ColorLimits',[0 0.7]);

EstimationMatrixPlotting(FC1,FC2)


%% Calculate Autocorrelation
sec = 3;
acf = zeros(N,sec/dt);
acf_sym = zeros(N,sec/dt);
for i = 1:N
    acf(i,:) = autocorr(S(i, 10/dt:end), 'NumLags', sec/dt-1);
    acf_sym(i,:) = autocorr(S2(i, 10/dt:end), 'NumLags', sec/dt-1);
end
%% Calculate Timescale
timeConstants = zeros(N, 1);
timeConstants_sym = zeros(N, 1);
for i = 1:N
        timeConstants(i) = calculate_time_constant(acf(i, :), dt);
        timeConstants_sym(i) = calculate_time_constant(acf_sym(i, :), dt); 
end

figure;
hold on;
bar(1:N, timeConstants, 'FaceColor',[33,49,80]./256, 'EdgeColor',[33,49,80]./256);
bar(1:N, timeConstants_sym, 'FaceColor', [199,35,54]./256, 'EdgeColor', [199,35,54]./256);
xticks(1:1:N);
    xticklabels({'V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc',...
        '7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2',...
        'STPi','ProM','F7','8B','STPr','24c'});
ylabel('Timescale');
alpha(0.7)
hold off;

function tau = calculate_time_constant(acf, dt)
    % Calculate timescale from autocorrelation
    startIndex = find(acf >= 0.9, 1, 'first');
    endIndex = find(acf <= 0.1, 1, 'first');
    if isempty(startIndex) || isempty(endIndex) || startIndex >= endIndex
        tau = NaN;
        return;
    end
    valid_x = (startIndex:endIndex) * dt;
    valid_acf = acf(startIndex:endIndex);
    
    ft = fittype('a * exp(-x/tau)', 'independent', 'x', 'dependent', 'y');
    opts = fitoptions(ft);
    opts.StartPoint = [max(valid_acf), 10]; 
    [fitresult, ~] = fit(valid_x', valid_acf', ft, opts);
    
    tau = fitresult.tau; 
end
