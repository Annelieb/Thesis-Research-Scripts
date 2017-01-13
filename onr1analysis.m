%% onr1
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
load('Subject_Accuracy_onr1.mat')
%% Subject_Accuracies onr1
plot(10:10:60, Subject_Accuracy)
legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12')
xlabel('Trial Length [s]')
ylabel('Accuracy [%]')
title('Subject Accuracies')

%% Mean & std
SAM = squeeze(mean(Subject_Accuracy,1));
SAS = squeeze(std(Subject_Accuracy,1));
figure
plot(10:10:60, SAM, 'k')
hold on
plot(10:10:60, SAM+SAS, '--k')
plot(10:10:60, SAM-SAS, '--k')
xlabel('Trial Length [s]')
ylabel('Accuracy [%]')
title('Averaged Accuracy over 12 Subjects')
axis([10 60 75 100.3])
hold off

%% Prior plots together
figure
plot(10:10:60, Subject_Accuracy)
hold on
plot(10:10:60, SAM, 'k', 'LineWidth', 2)
plot(10:10:60, SAM+SAS, '--k', 'LineWidth', 2)
plot(10:10:60, SAM-SAS, '--k', 'LineWidth', 2)
xlabel('Trial Length [s]')
ylabel('Accuracy [%]')
title('Averaged Accuracy over 12 Subjects')
axis([10 60 75 100.3])
legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'mean')
hold off

%% Compare with the pilot studies
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Investigate LP and HP values')
load('Result_LP1_HP8.mat')
Acc_pilot = squeeze(Subject_Accuracy(:,8,3,:));
PM = mean(Acc_pilot,2);
sPM = std(Acc_pilot,[],2);
PMn2 = mean(Acc_pilot(:,[1 3 4]),2);
sPMn2 = std(Acc_pilot(:,[1 3 4]),[],2);


figure
plot(10:10:60, SAM, 'k')
hold on
plot(10:10:60, PM, 'b')
plot(10:10:60, PMn2, 'r')
%plot(10:10:60, PM + sPM, '--b')
%plot(10:10:60, PM - sPM, '--b')
%plot(10:10:60, PMn2 + sPMn2, '--r')
%plot(10:10:60, PMn2 - sPMn2, '--r')
plot(10:10:60, SAM+SAS, '--k')
plot(10:10:60, SAM-SAS, '--k')
xlabel('Trial Length [s]')
ylabel('Accuracy [%]')
title('Compare Final Experiment with Pilot Experiment')
legend('Final Experiment', 'Pilot Experiment', 'Pilot Experiment without P2') 
axis([10 60 70 100.3])
set(gca, 'Xtick', [10:10:60])
hold off

% ttest
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
load('Subject_Accuracy_onr1.mat')
for trl = 1:6
    [h(trl) p(trl)] = ttest2(Subject_Accuracy(trl,:), Acc_pilot(trl,:), 'Tail', 'right')
    meandiff = mean(SAM - PM')
    sdiff = std(SAM - PM')
    
    [hn2(trl) pn2(trl)] = ttest2(Subject_Accuracy(trl,:), Acc_pilot(trl,[1 3 4]), 'Tail', 'right')
    meandiffn2 = mean(SAM - PMn2')
    sdiffn2 = std(SAM - PMn2')
end
%% onr1 difference between 'best_decoder' & 'leave one out'

dFBA = squeeze(FB_Accuracy(:,:,2)-FB_Accuracy(:,:,1))
MdFBA = squeeze(mean(dFBA,1))
