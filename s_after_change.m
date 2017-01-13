%% First second of 'changes' 
subjects = 1:12;
sort_trials = cell(12, 4, 4, 4);

for S = subjects
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Thresholds_Subject' num2str(S) '.mat'])
    for session = 1:4
        if S == 3 && session == 1
            load(['FB_S3_ses1_R_corrected_diff&Acc'], 'difference')
        end
        try
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_right.mat'], 'difference');
        catch
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_left.mat'], 'difference')
        end
        if length(difference) > 39
            difference(40:end) = []
        end
        for trial = 1:length(difference)
            if difference(trial) > th1;
                kleur{S,session}(trial) = 4;
            elseif difference(trial) > th2;
                kleur{S,session}(trial) = 3;
            elseif difference(trial) > th3;
                kleur{S,session}(trial) = 2;
            else;
                kleur{S,session}(trial) = 1;
            end;
        end
        for trial = 2:length(difference)
            k = kleur{S,session}(trial-1);
            l = kleur{S,session}(trial);
            sort_trials{S, session, k,l} = [sort_trials{S, session , k,l}, trial]
        end 
    end
end
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save('Sorted_trials', 'sort_trials')

%% Get the EEG data 100 previous samples and 200 after sections... 
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load('Sorted_trials')

for S = 1:12
    num = 0;
    for session = 1:4
        cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
        try
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_right.mat'], 'data2');
        catch
            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_left.mat'], 'data2')
        end
        data2(cellfun('isempty',data2)) = []; %Remove empty cells
        if size(data2,2)>39
            data2(40:end) = [];
        end
        for k = 1:4
            for l = 1:4
                %num=0;
                noftrials = length(sort_trials{S,session,k,l})
                for trl = 1: noftrials
                    num = num+1;
                    trial = sort_trials{S,session,k,l}(trl)+1; % +1 Because we want the 10s AFTER the change was shown to the subject!
                    if trial < 40
                        EEG_sorted_12{S, k,l, num} = [data2{trial-1} data2{trial}];
                        z = size(data2{trial-1},2)
                        for ch = 1:24
                            EEG_sorted_filtered{S, k,l,num}(ch,:) =  Rbp(1, 20, 6, 500, squeeze(EEG_sorted_12{S, k,l,num}(ch,:)));
                        end
                        EEG_sorted_cut{S,k,l}(:,:,num) = EEG_sorted_filtered{S,k,l,num}(:,z-100 : z + 200)
                        earref = (EEG_sorted_cut{S,k,l}(18,:,num)+ EEG_sorted_cut{S,k,l}(19,:,num))/2
                        for ch = 1:24
                            EEG_sorted_reref {S,k,l}(ch,:,num)= EEG_sorted_cut{S,k,l}(ch,:,num)-earref
                            first100 = mean(squeeze(EEG_sorted_reref {S,k,l}(ch,1:100,num)))
                            EEG_sorted_final{S,k,l}(ch,:,num) = squeeze(EEG_sorted_reref {S,k,l}(ch,:,num)) - first100
                        end
                    end
                end
                
            end
        end
    end
end

EEG_sorted_FB = EEG_sorted_final([1 2 5 6 9 10],:,:);
EEG_sorted_NFB = EEG_sorted_final([3 4 7 8 11 12],:,:);
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save(['EEG_sorted'], 'EEG_sorted_12', 'EEG_sorted_reref', 'EEG_sorted_final', 'EEG_sorted_filtered', 'EEG_sorted_cut', 'EEG_sorted_FB', 'EEG_sorted_NFB', '-v7.3')

%% Put all together for comparing FB and NFB
%cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
%load(['EEG_sorted'])

EEG_sorted_FB_all = cell(4,4);
EEG_sorted_NFB_all = cell(4,4);
EEG_sorted_FB_rtg = []; EEG_sorted_FB_gtr = []; % red to green and green to red
EEG_sorted_NFB_rtg = []; EEG_sorted_NFB_gtr = [];
EEG_sorted_FB_nochange = []; EEG_sorted_NFB_nochange = [];

for k = 1:4
    for l = 1:4
        for S = 1:6
            if isempty(EEG_sorted_FB{S,k,l})
            else
                EEG_sorted_FB_all{k,l} = cat(3,EEG_sorted_FB_all{k,l}, squeeze(EEG_sorted_FB{S,k,l}));
                if k < 3 && l > 2
                    EEG_sorted_FB_rtg = cat(3,EEG_sorted_FB_rtg, squeeze(EEG_sorted_FB{S,k,l}));
                elseif k > 2 && l < 3
                    EEG_sorted_FB_gtr = cat(3,EEG_sorted_FB_gtr, squeeze(EEG_sorted_FB{S,k,l}));
                elseif k == l
                    EEG_sorted_FB_nochange = cat(3,EEG_sorted_FB_nochange, squeeze(EEG_sorted_FB{S,k,l}));
                end
            end
            if isempty(EEG_sorted_NFB{S,k,l})
            else
                EEG_sorted_NFB_all{k,l} = cat(3,EEG_sorted_NFB_all{k,l}, squeeze(EEG_sorted_NFB{S,k,l}));
                if k < 3 && l > 2
                    EEG_sorted_NFB_rtg = cat(3,EEG_sorted_NFB_rtg, squeeze(EEG_sorted_NFB{S,k,l}));
                elseif k > 2 && l < 3
                    EEG_sorted_NFB_gtr = cat(3,EEG_sorted_NFB_gtr, EEG_sorted_NFB{S,k,l});
                elseif k == l
                    EEG_sorted_NFB_nochange = cat(3,EEG_sorted_NFB_nochange, squeeze(EEG_sorted_NFB{S,k,l}));
                end
            end
        end
    end
end

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save(['More_sorted_EEG.mat'], 'EEG_sorted_FB_all', 'EEG_sorted_FB_rtg', 'EEG_sorted_FB_gtr','EEG_sorted_FB_nochange','EEG_sorted_NFB_all', 'EEG_sorted_NFB_rtg', 'EEG_sorted_NFB_gtr', 'EEG_sorted_NFB_nochange', 'EEG_sorted_FB', 'EEG_sorted_NFB')

%% Averaging & plots

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load(['More_sorted_EEG.mat'], 'EEG_sorted_FB_all', 'EEG_sorted_FB_rtg', 'EEG_sorted_FB_gtr','EEG_sorted_FB_nochange','EEG_sorted_NFB_all', 'EEG_sorted_NFB_rtg', 'EEG_sorted_NFB_gtr', 'EEG_sorted_NFB_nochange', 'EEG_sorted', 'EEG_sorted_FB', 'EEG_sorted_NFB')


    keep = sum(squeeze(EEG_sorted_FB_rtg(1,:,:)));
    i = find(keep);
        EEG_sorted_FB_rtg(:,:,i) = [];
FB_rtg = squeeze(mean(EEG_sorted_FB_rtg,3));
st_FB_rtg = std(EEG_sorted_FB_rtg,[],3);
FB_gtr = squeeze(mean(EEG_sorted_FB_gtr,3));
st_FB_gtr = std(EEG_sorted_FB_gtr,[],3);
FB_nochange = squeeze(mean(EEG_sorted_FB_nochange,3));
st_FB_nochange = std(EEG_sorted_FB_nochange,[],3);

% NFB_rtg = squeeze(mean(EEG_sorted_NFB_rtg,3));
% st_NFB_rtg = std(EEG_sorted_NFB_rtg,[],3);
% NFB_gtr = squeeze(mean(EEG_sorted_NFB_gtr,3));
% st_NFB_gtr = std(EEG_sorted_NFB_gtr,[],3);
% NFB_nochange = squeeze(mean(EEG_sorted_NFB_nochange,3));
% st_NFB_nochange = std(EEG_sorted_NFB_nochange,[],3);
% 
for ch = 1:24
    figure
    plot(-99:200, FB_rtg(ch,:), 'g')
    hold on
    plot(-99:200, FB_gtr(ch,:), 'r')
    plot(-99:200, FB_nochange(ch,:), 'k')
    title(['Channel ' num2str(ch)])
    hold off
end