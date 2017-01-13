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

figure
plot(10:10:60, SAM, 'k')
hold on
plot(10:10:60, PM, 'b')
plot(10:10:60, PM + sPM, '--b')
plot(10:10:60, PM - sPM, '--b')
plot(10:10:60, SAM+SAS, '--k')
plot(10:10:60, SAM-SAS, '--k')
xlabel('Trial Length [s]')
ylabel('Accuracy [%]')
title('Compare Final Experiment with Pilot Experiment')
legend('Final Experiment', 'Pilot Experiment') 
%axis([10 60 75 100.3])
hold off

%% onr1 difference between 'best_decoder' & 'leave one out'

dFBA = squeeze(FB_Accuracy(:,:,2)-FB_Accuracy(:,:,1))
MdFBA = squeeze(mean(dFBA,1))
