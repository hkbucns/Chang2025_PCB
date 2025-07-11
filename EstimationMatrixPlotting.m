% This function is build for visualizing the distance and correlation
% between Jacobian matrix A and B.

% Input:
% A: Objective Variable,
% B: Estimated Variable.

% Output:
% SSE:          Normalized Sum-Squared Error,
% Corr:         Correlation of elements between non-diagonal A and B,
% Corr_nonzero: Correlation of non-zero and non-diagonal elements,
% A_vec_non:    Non-zero and non-diagonal vectorlized elements of A,
% B_vec_non:    Non-zero and non-diagonal vectorlized elements of B.

function [SSE,Corr,Corr_nonzero,A_vec_non,B_vec_non] = EstimationMatrixPlotting(A,B)

N = length(A);

%% Calculating Normalized Sum-Squared Error
SSE = norm(A-B)/norm(A); 

%% Plotting Elements Comparison
% Annotation: Diagonal elements are removed to reduce the effect of large value on correlation.
%             Hidden is the ttest between estimated and analytical Jacobian.

A_vec = reshape(A-diag(diag(A)),N^2,1);
B_vec = reshape(B-diag(diag(B)),N^2,1); 
Corr = corr(A_vec,B_vec);

% [Ttest,p] = ttest(A_vec,B_vec,"Alpha",0.001)

figure(3)
scatter(A_vec,B_vec,'MarkerEdgeColor',[243/256,169/256,147/256],...
        'MarkerFaceColor',[243/256,169/256,147/256]);
hold on
plot([0,max(A_vec)],[0,max(A_vec)],'--','Color','k');
hold off
set(gca,'box','off');
% title('Estimated SC against Analytical SC');
xlabel('Ground Truth SC');
ylabel('Estimated SC');
% xlabel('w-I Jacobian');
% ylabel('\tau-b Jacobian');
xticks([0 0.8]);
yticks([0 0.8]);
alpha(0.8);
saveas(gcf,'SCPerformance.png');

%% Plotting Non-zero Elements Comparison
% Annotation: Non-zero elements Reconstruction Performance.
%             Hidden is the ttest between non-zero estimated and analytical Jacobian.

Index = (A_vec~=0);
A_vec_non = A_vec(Index);
B_vec_non = B_vec(Index);
Corr_nonzero = corr(A_vec_non,B_vec_non);
% SSE = norm(A_vec_non-B_vec_non)/norm(A_vec_non);
% [Ttest_non,p_non] = ttest(A_vec_non,B_vec_non)

figure(4)
scatter(A_vec_non,B_vec_non,'MarkerEdgeColor',[243/256,169/256,147/256],...
        'MarkerFaceColor',[243/256,169/256,147/256]);
hold on
plot([0,max(A_vec)],[0,max(A_vec)],'--','Color','k');
hold off
set(gca,'box','off');
% title('Estimated SC against Analytical SC (Non zero)');
xlabel('Ground Truth SC');
ylabel('Estimated SC');
xticks([0 0.8]);
yticks([0 0.8]);
alpha(0.8);
saveas(gcf,'NonZeroSCPerformance.png');

%% Plotting single direction coupling

%% Plotting Degree Distributions (In-Out)
% Annotation: Unfinished.

% A_Index = ((A-diag(diag(A)))~=0);
% B_Index = ((B-diag(diag(B)))~=0);
% A_Degree_In = sum(A_Index,2);
% A_Degree_Out = sum(A_Index,1);
% B_Degree_In = sum(B_Index,2);
% B_Degree_Out = sum(B_Index,1);
% 
% [A_Degree_In, A_In_Index] = sort(A_Degree_In);
% A_Degree_Out = A_Degree_Out(A_In_Index);
% 
% 
% figure(3)

%% 


end