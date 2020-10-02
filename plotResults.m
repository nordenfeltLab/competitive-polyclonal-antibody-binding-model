function plotResults(TotalBinding, FabBinding, FcBinding, conc)
%% Plotting Binding Probability Curve
figure(2)
hold on

%plot(conc, FabBinding, conc, FcBinding)
%plot(conc, FcBinding)
plot(conc, FabBinding, conc, FcBinding, 'LineWidth', 1)
%plot(cone,totexp,'*')
%plot(conc,FabBinding, 'LineWidth', 2)
%set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([10 10000]);
ylim([0 1.5]);
title('IvIgG Sf370');
xlabel('Total IgG concentration')
ylabel('Bound IgG')
%legend('Non competitive Fab binding','Competitive Fab binding','Non competitive Fc binding','Competitive Fc binding')
legend('Theoretical Fab-binding','Theoretical Fc-binding')
%legend('Fab binding','Fc binding')
%legend('Competitive Fab')
%legend('Fab binding experimental', 'Fc binding experimental', 'Total binding experimental','Fab binding computed', 'Fc binding computed', 'Total binding computed')
%legend('Fab binding experimental', 'Fab binding computed')
