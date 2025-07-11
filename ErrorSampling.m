%% Parameters setting

load('SC.mat');
SC = fln;
N = length(SC);
T = 50000; %s, ;
dt = 0.01; %10ms, ;
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


%% Repeatition
repeat_times = 10;
Corr_sampling_exp = zeros(repeat_times,11);
Corr_sampling_J = zeros(repeat_times,11);

for repeat = 1:repeat_times
%% Reduced Wong-Wang Model Simulation

[S,eta] = dMFM(SC, dt, T, w, I, G, sigma);
S_star = mean(S,2);
x_star = w.*J.*S_star+G.*J.*SC*S_star+I;

SC_sym = (SC + SC')/2;
dH_val = dH(x_star); % reveal the value of dH/dt at x_star


%% Jacobian Matrix Calculation

Jacob = zeros(N);
y_ana = zeros(N,1);
for i = 1:N
    y_ana(i) = 1./(gamma*G*J*(1-S_star(i))*dH(x_star(i)));
    for j = 1:N
        if i == j
            Jacob(i,j) = -1/(tau*(1-S_star(i)))+w(i)*gamma*J*(1-S_star(i))*dH(x_star(i));
        else
            Jacob(i,j) = gamma*G*J*(1-S_star(i))*SC(i,j)*dH(x_star(i));
        end
    end
end
clear i j 

x = S;

[N,M]=size(x);
L = zeros(N);
dx = (x(:,2:5e4)-x(:,1:5e4-1))./dt;
C2 = cov([dx' x(:,1:5e4-1)']); % Calculating Partial Cov
C1 = C2(1:N,N+1:2*N);
C = cov(x(:,1:5e4-1)');
L(:,:) = C1/C;

[y_st,w_recon,C_recon] = RevealHHetero1(SC_sym,L,S_star,tau,gamma,G, J);
C_recon = C_recon-diag(diag(C_recon));
dH_st = 1./(gamma*G*J.*(1-S_star).*y_st);

A_vec = reshape(Jacob-diag(diag(Jacob)),N^2,1);
B_vec = reshape(L-diag(diag(L)),N^2,1);
Corr_sampling_exp(repeat,1) = corr(A_vec,B_vec);

control = A_vec;
Corr_sampling_J(repeat,1) = corr(control,B_vec);


color1 = [33,49,80]./256;
color2 = [199, 35, 54] ./ 256;

figure(1) % Plotting the Figure 3A
subplot(2,6,1)
scatter(A_vec,B_vec,'MarkerEdgeColor',color1,...
        'MarkerFaceColor',color1);
hold on
plot([0,max(A_vec)],[0,max(A_vec)],'--','Color','k','LineWidth',1.5);
hold off
set(gca,'box','off');
xticks([]);
yticks([]);
ylim([-0.05 max(B_vec)]);
alpha(0.5);
xlabel('0.01s (Unsampled)');
set(gca, 'FontName', 'Arial')

for i=2:11
    rate = 10*(i-1);
    x = downsample(S',rate)';
    dT = dt*rate;
    [N,M]=size(x);
    L = zeros(N);
    M_sample = 5e4;
    dx = (x(:,2:M_sample)-x(:,1:M_sample-1))./dT;
    C2 = cov([dx' x(:,1:M_sample-1)']); % Calculating Partial Cov
    C1 = C2(1:N,N+1:2*N);
    C = cov(x(:,1:M_sample-1)');
    L(:,:) = C1/C;
    Jacob_sampling_exp = (expm(dT*Jacob)-eye(N))./dT;
    J_est = logm(L.*dT+eye(N))./dT;

    B_vec = reshape(L-diag(diag(L)),N^2,1);
    C_vec = reshape(Jacob_sampling_exp-diag(diag(Jacob_sampling_exp)),N^2,1);
    J_vec = reshape(J_est-diag(diag(J_est)),N^2,1);

    Corr_sampling_exp(repeat,i) = corr(C_vec,B_vec);
    Corr_sampling_J(repeat,i) = corr(control,J_vec);
    subplot(2,6,i)
    scatter(C_vec,B_vec,'MarkerEdgeColor',color1,...
            'MarkerFaceColor',color1);
    hold on
    plot([0,max(C_vec)],[0,max(C_vec)],'--','Color','k','LineWidth',1.5);
    hold off
    set(gca,'box','off');
    xticks([]);
    yticks([]);
    ylim([-0.05 max(B_vec)]);
    alpha(0.5);
    xlabel(sprintf('0.%i s',i-1));
    set(gca, 'FontName', 'Arial')

end
end

%% Plotting the Figure 3B

figure(2)
plot(0:10:100,mean(Corr_sampling_exp,1),'Color',color1, 'LineWidth', 2);
hold on
plot(0:10:100,mean(Corr_sampling_J,1),'Color',color2, 'LineWidth', 2);hold on
patch([0:10:100 fliplr(0:10:100)],...
    [(mean(Corr_sampling_J,1)-std(Corr_sampling_J,0,1)) fliplr((mean(Corr_sampling_J,1)+std(Corr_sampling_J,0,1)))],...
    color2,'edgealpha', '0', 'facealpha', '.2')
patch([0:10:100 fliplr(0:10:100)],...
    [(mean(Corr_sampling_exp,1)-std(Corr_sampling_exp,0,1)) fliplr((mean(Corr_sampling_exp,1)+std(Corr_sampling_exp,0,1)))],...
    color1,'edgealpha', '0', 'facealpha', '.2')
hold off
legend('exponential J','J');
set(gca,'box','off');
xticks(0:10:100);
yticks([0.8 1]);
ylim([0 1]);
xlabel('Sampling Steps');
ylabel('${\rm Corr({J},L)}$','interpreter','latex','FontSize',14);
set(gca, 'FontName', 'Arial')

