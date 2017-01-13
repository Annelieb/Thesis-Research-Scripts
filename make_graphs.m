%% ANALYSIS & GRAPHS OF RESULTS

cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\Results 1, 2, 3, 4')
for sessions = 1:4
    load(['Results_' num2str(sessions) '.mat'])
    eval(['Accuracy' num2str(sessions) '= Accuracy;'])
    eval(['CA' num2str(sessions) '= CA;'])
    eval(['CUA' num2str(sessions) '= CUA;'])
    Accuracy = []; CUA = []; CA = [];
end

for sessions = 1:4
    load(['Results_per_type' num2str(sessions) '.mat'])
    eval(['Accuracy_per_type' num2str(sessions) '= Accuracy_per_type;'])
    eval(['decoders' num2str(sessions) '= decoders;'])
    types(sessions) = noftypes;
    Accuracy_per_type = []; decoders = []; noftypes =[];
end

% Intermediate500
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate')
load('Intermediate500')

%% Calculate some needed matrices
for sessions = 1:4
    eval(['Sub_mean_' num2str(sessions) '= squeeze(mean(Accuracy' num2str(sessions) ', 3));'])
    eval(['Sub_std_' num2str(sessions) '= std(Accuracy' num2str(sessions) ',0,3);'])
    for type = 1:types(sessions)
        eval(['Sub_mean_pt_' num2str(sessions) '{type}= squeeze(mean(Accuracy_per_type' num2str(sessions) '{type}, 3));'])
        eval(['Sub_std_pt_' num2str(sessions) '{type}= std(Accuracy_per_type' num2str(sessions) '{type},0,3);'])
    end
end

%%
All_S_Acc_1 = []; All_S_Acc_2 = [];All_S_Acc_3 = []; All_S_Acc_4 = [];
All_S_Acc_pt_1 = {} ; All_S_Acc_pt_2 = {};All_S_Acc_pt_3 = {}; All_S_Acc_pt_4 = {};

for S = 1:4
    for sessions = 1:4
    eval(['All_S_Acc_' num2str(sessions) '= cat(2,All_S_Acc_' num2str(sessions) ',squeeze(Accuracy' num2str(sessions) '(:,S,:)));'])
%         for type = 1:types(sessions)
%             eval(['All_S_Acc_pt_' num2str(sessions) '{type}= cat(2,All_S_Acc_pt_' num2str(sessions) '{type},squeeze(Accuracy_per_type' num2str(sessions) '{type}(:,S,:)));'])
%         end
    end
end

%% Because above doesn't work
for sessions = 1:4
    for type = 1:types(sessions)
        eval(['Tot_mean_pt_' num2str(sessions) '{type}= squeeze(mean(Sub_mean_pt_' num2str(sessions) '{type}, 2));'])
        eval(['Tot_std_pt_' num2str(sessions) '{type}= std(Sub_mean_pt_' num2str(sessions) '{type},0,2);'])
    end
end

%%
for sessions = 1:4
    eval(['Tot_mean_' num2str(sessions) '= squeeze(mean(All_S_Acc_' num2str(sessions) ', 2));'])
    eval(['Tot_std_' num2str(sessions) ' = std(All_S_Acc_' num2str(sessions) ',0,2);'])
%     for type = 1:types(sessions)
%         eval(['Tot_mean_pt_' num2str(sessions) '{type}= squeeze(mean(All_S_Acc_pt_' num2str(sessions) '{type}, 2));'])
%         eval(['Tot_std_pt_' num2str(sessions) '{type} = std(All_S_Acc_pt_' num2str(sessions) '{type},0,2);'])
%     end
end

Tot_mean_I = squeeze(mean(Intermediate500, 2));
Tot_std_I = std(Intermediate500,0,2);

%%
% %% Not 2
% All_S_Acc_1 = []; All_S_Acc_2 = [];All_S_Acc_3 = []; All_S_Acc_4 = [];
% for S = [1, 3, 4]
%     All_S_Acc_1 = cat(2,All_S_Acc_1,squeeze(Accuracy1(:,S,:)));
%     All_S_Acc_2 = cat(2,All_S_Acc_2,squeeze(Accuracy2(:,S,:)));
%     All_S_Acc_3 = cat(2,All_S_Acc_3,squeeze(Accuracy3(:,S,:)));
%     All_S_Acc_4 = cat(2,All_S_Acc_4,squeeze(Accuracy4(:,S,:))); 
% end
% 
% Tot_mean_1 = squeeze(mean(All_S_Acc_1, 2));
% Tot_std_1 = std(All_S_Acc_1,0,2);
% 
% Tot_mean_2 = squeeze(mean(All_S_Acc_2, 2));
% Tot_std_2 = std(All_S_Acc_2,0,2);
% 
% Tot_mean_3 = squeeze(mean(All_S_Acc_3, 2));
% Tot_std_3 = std(All_S_Acc_3,0,2);
% 
% Tot_mean_4 = squeeze(mean(All_S_Acc_4, 2));
% Tot_std_4 = std(All_S_Acc_4,0,2);
% 
% Tot_mean_I = squeeze(mean(Intermediate500(:, [1, 3, 4]), 2));
% Tot_std_I = std(Intermediate500(:, [1, 3, 4]),0,2);

 %% Accuracies per subject, averaged over decoders
% 
% figure
% 
% plot(10:10:60, Sub_mean_2(:,1), 'r')
% hold on
% plot(10:10:60, Sub_mean_2(:,2), 'b')
% plot(10:10:60, Sub_mean_2(:,3), 'g')
% plot(10:10:60, Sub_mean_2(:,4), 'm')
% 
% plot(10:10:60, Sub_mean_2(:,1) + Sub_std_2(:,1), 'r--')
% plot(10:10:60, Sub_mean_2(:,1) - Sub_std_2(:,1), 'r--')
% 
% plot(10:10:60, Sub_mean_2(:,2) + Sub_std_2(:,1), 'b--')
% plot(10:10:60, Sub_mean_2(:,2) - Sub_std_2(:,1), 'b--')
% 
% plot(10:10:60, Sub_mean_2(:,3) + Sub_std_2(:,1), 'g--')
% plot(10:10:60, Sub_mean_2(:,3) - Sub_std_2(:,1), 'g--')
% 
% plot(10:10:60, Sub_mean_2(:,4) + Sub_std_2(:,1), 'm--')
% plot(10:10:60, Sub_mean_2(:,4) - Sub_std_2(:,1), 'm--')
% 
% % Title etc.
% title('Averaged accuracies over decoders')
% legend('S1', 'S2' ,'S3', 'S4')
% xlabel('Trial Lengths [s]')
% ylabel('Accuracies')
% 
% hold off

%% Averages over decoders compared to other results per subject
for S = 1:4
    figure
    
    plot(10:10:60, Sub_mean_1(:,S), 'b')
    hold on
    
    plot(10:10:60, Sub_mean_2(:,S), 'r')
    plot(10:10:60, Sub_mean_3(:,S), 'g')
    plot(10:10:60, Sub_mean_4(:,S), 'm')
    plot(10:10:60, Intermediate500(:,S), 'k')
    
    plot(10:10:60, Sub_mean_1(:,S) + Sub_std_1(:,S), 'b--')
    plot(10:10:60, Sub_mean_1(:,S) - Sub_std_1(:,S), 'b--')
        
    plot(10:10:60, Sub_mean_2(:,S) + Sub_std_2(:,S), 'r--')
    plot(10:10:60, Sub_mean_2(:,S) - Sub_std_2(:,S), 'r--')
    
    plot(10:10:60, Sub_mean_3(:,S) + Sub_std_3(:,S), 'g--')
    plot(10:10:60, Sub_mean_3(:,S) - Sub_std_3(:,S), 'g--')
    
    plot(10:10:60, Sub_mean_4(:,S) + Sub_std_4(:,S), 'm--')
    plot(10:10:60, Sub_mean_4(:,S) - Sub_std_4(:,S), 'm--')
    
    % Title etc.
    title(['Compare averaged accuracies over decoderes: Subject ' num2str(S)])
    legend('1 Session', '2 Sessions', '3 Sessions', '4 Sessions', 'Leave one out')
    xlabel('Trial Lengths [s]')
    ylabel('Accuracies')

    hold off
end

%% Averaged per type for each subject and each number of sessions
for S = 1:4
    for sessions = 1:4
        figure
        for type = 1:types(sessions)
            eval(['plot(10:10:60, Sub_mean_pt_' num2str(sessions) '{type}(:,S))'])
            hold on
        end
    
    plot(10:10:60, Intermediate500(:,S), 'k')
    eval(['plot(10:10:60, Sub_mean_' num2str(sessions) '(:,S), ''k-*'')'])
    title(['Compare types: Subject ' num2str(S) ', Number of sessions = ' num2str(sessions)])
    xlabel('Trial Lengths [s]')
    ylabel('Accuracies')
    
    if sessions == 1
        legend('Type 1', 'Type 2', 'leave one out', 'Avg of all') 
    elseif sessions == 2
        legend('Type 1', 'Type 2', 'Type 3', 'leave one out', 'Avg of all')    
    elseif sessions == 3
        legend('Type 1', 'Type 2', 'Type 3', 'Type 4', 'leave one out', 'Avg of all')     
    elseif sessions == 4
        legend('Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5', 'leave one out', 'Avg of all') 
    end
    hold off
    end
end
    
%% Averaged per type over the subjects for each number of sessions

for sessions = 1:4

    figure
    for type = 1:types(sessions)
        eval(['plot(10:10:60, Tot_mean_pt_' num2str(sessions) '{type})'])
        hold on
    end
    
    plot(10:10:60, Tot_mean_I, 'k')
    eval(['plot(10:10:60, Tot_mean_' num2str(sessions) ', ''k-*'')'])
    title(['Averaged accuracies, compare types. Number of sessions = ' num2str(sessions)])
    xlabel('Trial Lengths [s]')
    ylabel('Accuracies')
    
    if sessions == 1
        legend('Type 1', 'Type 2', 'leave one out', 'Avg of all') 
    elseif sessions == 2
        legend('Type 1', 'Type 2', 'Type 3', 'leave one out', 'Avg of all')    
    elseif sessions == 3
        legend('Type 1', 'Type 2', 'Type 3', 'Type 4', 'leave one out', 'Avg of all')     
    elseif sessions == 4
        legend('Type 1', 'Type 2', 'Type 3', 'Type 4', 'Type 5', 'leave one out', 'Avg of all') 
    end
    hold off
end
%% Averaged over all

figure


    plot(10:10:60, Tot_mean_1, 'b')
    hold on
    plot(10:10:60, Tot_mean_2, 'r')
    plot(10:10:60, Tot_mean_3, 'g')
    plot(10:10:60, Tot_mean_4, 'm')
    plot(10:10:60, Tot_mean_I, 'k')
    
    plot(10:10:60, Tot_mean_1 + Tot_std_1, 'b--')
    plot(10:10:60, Tot_mean_1 - Tot_std_1, 'b--')
    
    plot(10:10:60, Tot_mean_2 + Tot_std_2, 'r--')
    plot(10:10:60, Tot_mean_2 - Tot_std_2, 'r--')

    plot(10:10:60, Tot_mean_3 + Tot_std_3, 'g--')
    plot(10:10:60, Tot_mean_3 - Tot_std_3, 'g--')
        
    plot(10:10:60, Tot_mean_4 + Tot_std_4, 'm--')
    plot(10:10:60, Tot_mean_4 - Tot_std_4, 'm--')

    plot(10:10:60, Tot_mean_I + Tot_std_I, 'k--')
    plot(10:10:60, Tot_mean_I - Tot_std_I, 'k--')
    
    % Title etc.
    title('Averaged accuracies over decoders and subjects')
    legend('1 Session', '2 Sessions', '3 Sessions', '4 Sessions', 'Leave one out')
    xlabel('Trial Lengths [s]')
    ylabel('Accuracies')

    hold off
    
    %% All decoders per subject
    for S = 1:4
    figure
    plot(10:10:60, squeeze(Accuracy2(:,S,:)))
    
    % Title etc.
    title(['Accuracies per decorder: Subject ' num2str(S)])

    xlabel('Trial Lengths [s]')
    ylabel('Accuracies')
    end
