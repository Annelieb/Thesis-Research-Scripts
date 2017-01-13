% GRAPHS OF DECODERS WITH DIFFERENT TRL
% Load the data
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\different trial lengths for evaluation and decoders')
for sessions = 1:4
    load(['Results_' num2str(sessions) '.mat'])
    eval(['Accuracy' num2str(sessions) '= Accuracy;'])
    eval(['CA' num2str(sessions) '= CA;'])
    eval(['CUA' num2str(sessions) '= CUA;'])
    Accuracy = []; CUA = []; CA = [];
end


%% Initial averageing
for sessions = 1:4
    eval(['mean_accuracy' num2str(sessions) ' = squeeze(mean(Accuracy' num2str(sessions) ',3));'])
    eval(['std_accuracy' num2str(sessions) ' = squeeze(std(Accuracy' num2str(sessions) ',0,3));'])
end

%% reorder the matrix
oldorder = [1 2 3 4 5 6 7 8];
neworder = [7 1 8 2 3 4 5 6];

for sessions = 1:4
    for  i= 1:8 % fix the first dimension
        eval(['mean_accuracy_n' num2str(sessions) '(' num2str(i) ', :, :) = mean_accuracy' num2str(sessions) '(' num2str(neworder(i)) ',:,:);'])
    end

    for i = 1:8 % fix the second dimension
        eval(['mean_accuracy_new' num2str(sessions) '(:,:,' num2str(i) ') = mean_accuracy_n' num2str(sessions) '(:,:,' num2str(neworder(i)) ');'])
    end
end
    

%% Show plots averaged over subjects
for sessions = 1:4
    eval(['subj_mean_accuracy' num2str(sessions) ' = squeeze(mean(mean_accuracy_new' num2str(sessions) ',2));'])
end
%%
for sessions = 1:4
    figure
    eval(['plot([5 10 15 20 30 40 50 60], subj_mean_accuracy' num2str(sessions) '(1:4, :)'')'])
    legend('Trial length 5', 'Trial length 10', 'Trial length 15', 'Trial length 20')
    title(['Decoders made from different trial lengths; Number of sessions used = ' num2str(sessions) ])
    xlabel('Trial length of decoder')
    ylabel('Accuracy')
end
%%
for sessions = 1:4
    eval(['plot([5 10 15 20 30 40 50 60], subj_mean_accuracy' num2str(sessions) '(1:4, :)'')'])
    hold on
end
legend('TRL5, ses1', 'TRL10, ses1','TRL15, ses1','TRL20, ses1', 'TRL5, ses2', 'TRL10, ses2','TRL15, ses2','TRL20, ses2','TRL5, ses3', 'TRL10, ses3','TRL15, ses3','TRL20, ses3','TRL5, ses4', 'TRL10, ses4','TRL15, ses4','TRL20, ses4')
    title(['Decoders made from different trial lengths'])
    xlabel('Trial length of decoder')
    ylabel('Accuracy')
    hold off
%% Plots averaged over subjects


%% Calculate differences between adjacent decoder lengths
for sessions = 1:4
    for f = 1:7
        eval(['diff' num2str(sessions) '(:,:,' num2str(f) ') = mean_accuracy_new' num2str(sessions) '(:,:,' num2str(f) ') - mean_accuracy_new' num2str(sessions) '(:,:,' num2str(f+1) ');'])
    end
end

%% Calculate differences with first
for sessions = 1:4
    for f = 1:8
        eval(['diff_from1_' num2str(sessions) '(:,:,f) = mean_accuracy_new' num2str(sessions) '(:,:,1) - mean_accuracy_new' num2str(sessions) '(:,:,' num2str(f) ');'])
    end
end