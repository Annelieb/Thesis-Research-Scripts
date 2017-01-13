%% Subject Accuracies samenvoegen
subjects = 1:12
for onr = 7
for S = subjects
    NSR = 20; HP = 8; LP = 1;
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
    load(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_onr' num2str(onr) '.mat'], 'subject_Accuracy')
    Subject_Accuracy(S, :) = subject_Accuracy
end

cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
save(['Subject_Accuracy_onr' num2str(onr) '.mat'], 'Subject_Accuracy')
end
%% Evaluate decoders: Feedback accuracy chosen decoder compared to leave one out
subjects = 1:12
trll = [10 20 30 40 50 60];
onr = 1; NSR = 20; LP = 1; HP = 8;
    
% load pilot
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load('Pilot_decoder_all.mat')
pilot_all = avg_dec;
load('Pilot_decoder_134.mat')
pilot_134 = avg_dec;


for S = subjects
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
    load(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_onr' num2str(onr) '.mat'], 'EE', 'Attended', 'Unattended', 'subject_difference', 'subject_CA', 'subject_CUA')
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Best_decoder_Subject' num2str(S) '.mat'], 'avg_dec')
    best = avg_dec;
    load(['All_decoders_Subject' num2str(S) '.mat'], 'decoders')
    four = decoders{4,1};
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
    load(['Decoder all other S' num2str(S) '.mat'], 'avg_dec')
    other = avg_dec;
    % load all other -> SCRIPT IS AF, MAAR MOET NOG BEREKEND WORDEN
    % other =
    
    for TRL = 1:length(trll)
        nofFB = size(subject_CA{TRL},2)/2;
        train = (nofFB +1 :nofFB*2);
        CAbest = []; CUAbest =[];
        CAfour = []; CUAfour =[];
        CApilot_all = []; CUApilot_all =[];
        CApilot_134 = []; CUApilot_134 =[];
        CAother = []; CUAother =[];

        for HH = 1: nofFB
            
            recenvbest = EE{TRL,train(end)} * best;
            recenvfour = EE{TRL,train(end)} * four;
            recenvpilot_all = EE{TRL,train(end)} * pilot_all{1,TRL};
            recenvpilot_134 = EE{TRL,train(end)} * pilot_134{1,TRL};
            recenvother = EE{TRL,train(end)} * other{1,TRL};


            CAbest(train(end)-nofFB) = corr2(recenvbest,Attended{TRL,train(end)}); % Attended
            CUAbest(train(end)-nofFB) = corr2(recenvbest,Unattended{TRL,train(end)}); % Unattended
            
            CAfour(train(end)-nofFB) = corr2(recenvfour,Attended{TRL,train(end)}); % Attended
            CUAfour(train(end)-nofFB) = corr2(recenvfour,Unattended{TRL,train(end)}); % Unattended
            
            CApilot_all(train(end)-nofFB) = corr2(recenvpilot_all,Attended{TRL,train(end)}); % Attended
            CUApilot_all(train(end)-nofFB) = corr2(recenvpilot_all,Unattended{TRL,train(end)}); % Unattended
 
            CApilot_134(train(end)-nofFB) = corr2(recenvpilot_134,Attended{TRL,train(end)}); % Attended
            CUApilot_134(train(end)-nofFB) = corr2(recenvpilot_134,Unattended{TRL,train(end)}); % Unattended
           
            CAother(train(end)-nofFB) = corr2(recenvother,Attended{TRL,train(end)}); % Attended
            CUAother(train(end)-nofFB) = corr2(recenvother,Unattended{TRL,train(end)}); % Unattended

            train = [train(end) train(1:end-1)];
        end
        % 'best decoder'
%         FB_CA{TRL}(S,1,:) = CA;
%         FB_CUA{TRL}(S,1,:) = CUA;
%         FB_difference{TRL}(S,1,:) = CA - CUA;
        FB_Accuracy(S,TRL,1) = (sum(CAbest-CUAbest> 0)/nofFB)*100
        
        % leave one out
        FB_CA{TRL}(S,2,:) = subject_CA{TRL}(nofFB +1 :nofFB*2);
        FB_CUA{TRL}(S,2,:) = subject_CUA{TRL}(nofFB +1 :nofFB*2);
        FB_difference{TRL}(S,2,:) = subject_difference{TRL}(nofFB +1 :nofFB*2);
        FB_Accuracy(S,TRL,2) = (sum(FB_difference{TRL}(S,2,:)> 0)/nofFB)*100
        
        % four
        FB_Accuracy(S,TRL,3) = (sum(CAfour-CUAfour> 0)/nofFB)*100
        
        % pilot all
        FB_Accuracy(S,TRL,4) = (sum(CApilot_all-CUApilot_all> 0)/nofFB)*100
        
        %pilot 134
        FB_Accuracy(S,TRL,5) = (sum(CApilot_134-CUApilot_134> 0)/nofFB)*100
        
        %other
        FB_Accuracy(S,TRL,6) = (sum(CAother-CUAother> 0)/nofFB)*100
    end
end

cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
save(['FB_compare_decoders.mat'], 'FB_Accuracy')

%% onr1 difference between 'best_decoder' & 'leave one out'
onr = 1
cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
load(['FB_Accuracies_leav1out_vs_bestdecoder.mat'], 'FB_CUA', 'FB_CA', 'FB_difference', 'FB_Accuracy')

dFBA = squeeze(FB_Accuracy(:,:,2)-FB_Accuracy(:,:,1))
MdFBA = squeeze(mean(dFBA,1))
SdFBA = std(dFBA, [], 1)
for TRL = 1:6
    [h,p] = ttest(FB_Accuracy(:,TRL,2), FB_Accuracy(:,TRL,1), 'Tail', 'right')
    P(TRL) = p
end

% Mean difference: plot + significance
mFBA = squeeze(mean(FB_Accuracy, 1))
sFBA = squeeze(std(FB_Accuracy, [], 1))


figure
plot(10:10:60, mFBA(:,1), 'b', 'LineWidth', 2)
hold on
plot(10:10:60, mFBA(:,2), 'k', 'LineWidth', 2)
plot(10:10:60, mFBA(:,1) + sFBA(:,1) , '--b')
plot(10:10:60, mFBA(:,1) - sFBA(:,1) , '--b')
plot(10:10:60, mFBA(:,2) + sFBA(:,2) , '--k')
plot(10:10:60, mFBA(:,2) - sFBA(:,2) , '--k')
hold off
title('Trained decoder vs Leave one out validation')
legend('Trained decoder', 'Leave one out')
xlabel('Trial Length [s]')
ylabel('Accuracy [%]')

% average difference
% dmFBA = mFBA(:,2)-mFBA(:,1)
% mdmFBA = mean(dmFBA)
% sdmFBA = std(dmFBA)
% 
% for TRL = 1:6
%     [h,p] = ttest(FB_Accuracy(:,TRL,2), FB_Accuracy(:,TRL,1), 'Tail', 'right')
%     P(TRL) = p
% end

%% Was the best decoder selected? 
onr = 1; NSR = 20; LP = 1; HP = 8;
trll = [10 20 30 40 50 60];

Accuracy_overall = cell(12,6);
Accuracy_T = cell(12,6);
Accuracy_T = cell(12,6);

for S = 1:12
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['All_decoders_Subject' num2str(S) '.mat'], 'decoders')
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
    load(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_onr' num2str(onr) '.mat'], 'EE', 'Attended', 'Unattended', 'subject_difference', 'subject_CA', 'subject_CUA')
    for dec = 1:4*6
        avg_dec = decoders{dec}
        if isempty(avg_dec)
        else
            for TRL = 1:length(trll)
                if isempty(Accuracy_overall{S, TRL})
                    Accuracy_overall{S,TRL} = zeros(4,6);
                    Accuracy_T{S,TRL} = zeros(4,6);
                    Accuracy_F{S,TRL} = zeros(4,6);
                end
                nofFB = size(subject_CA{TRL},2);
                train = (1:nofFB);
                CA = []; CUA =[];
                for HH = 1: nofFB

                    recenv = EE{TRL,train(end)} * avg_dec;

                    CA(train(end)) = corr2(recenv,Attended{TRL,train(end)}); % Attended
                    CUA(train(end)) = corr2(recenv,Unattended{TRL,train(end)}); % Unattended

                    train = [train(end) train(1:end-1)];
                end
                
                difference = CA - CUA;
                Accuracy_overall{S,TRL}(dec) = (sum(CA-CUA> 0)/nofFB)*100
                Accuracy_T{S,TRL}(dec)= (sum(CA(1:nofFB/2)-CUA(1:nofFB/2)> 0)/(nofFB/2))*100
                Accuracy_F{S,TRL}(dec)= (sum(CA(nofFB/2+1:end)-CUA(nofFB/2+1:end)> 0)/(nofFB/2))*100
            end
        end
    end
end

% Find the best decoders
for S = 1:12
    for TRL = 1:6
        [Max_Acc_overall{S,TRL}, I_all{S,TRL}] = max(Accuracy_overall{S, TRL}(:)); 
        [Max_Acc_T{S,TRL}, I_T{S,TRL}] = max(Accuracy_T{S,TRL}(:)); 
        [Max_Acc_F{S,TRL}, I_F{S,TRL}] = max(Accuracy_F{S,TRL}(:)); 
    end
end

%% Difference between before and after FB & NFB
subjects = 1:12
trll = [10 20 30 40 50 60];
onr = 1; NSR = 20; LP = 1; HP = 8;

for S = subjects
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
    load(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_onr' num2str(onr) '.mat'], 'subject_difference', 'subject_CA', 'subject_CUA')
    for TRL = 1:length(trll)
        helft = size(subject_CA{TRL},2)/2;
        
        % Training
        TF_CA{TRL}(S,1,:) = subject_CA{TRL}(1:helft);
        TF_CUA{TRL}(S,1,:) = subject_CUA{TRL}(1:helft);
        TF_difference{TRL}(S,1,:) = subject_difference{TRL}(1:helft);
        TF_Accuracy(S,TRL,1) = (sum(TF_difference{TRL}(S,1,:)>0)/helft)*100
        
        % Feedback
        TF_CA{TRL}(S,2,:) = subject_CA{TRL}(helft+1:end);
        TF_CUA{TRL}(S,2,:) = subject_CUA{TRL}(helft+1:end);
        TF_difference{TRL}(S,2,:) = subject_difference{TRL}(helft+1:end);
        TF_Accuracy(S,TRL,2) = (sum(TF_difference{TRL}(S,2,:)>0)/helft)*100
    end
end

cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
save(['Accuracies_TvsF.mat'], 'TF_CUA', 'TF_CA', 'TF_difference', 'TF_Accuracy')

%% Difference between T&F
% Above
onr = 1
cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
load(['Accuracies_TvsF.mat'])

% OR!!! Online
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
load('Online experimentes accuracies put together.mat', 'Session_Accuracy')
TF_Accuracy(:,1,1) = squeeze(mean(Session_Accuracy(:,1:4),2)); 
TF_Accuracy(:,1,2) = squeeze(mean(Session_Accuracy(:,5:8),2)); 


dTFA = squeeze(TF_Accuracy(:,:,1)-TF_Accuracy(:,:,2))
avgdTFA = mean(dTFA,1)
% paired ttest
    for T = 1: size(TF_Accuracy,2)
        [h p] = ttest(squeeze(TF_Accuracy(:,T,1)), squeeze(TF_Accuracy(:,T,2)), 'Alpha',0.05)
        result(T,:) = [h p]
    end

% split in feedback vs no feedback
NFBgroup = TF_Accuracy([3 4 7 8 11 12],:,:)
FBgroup = TF_Accuracy([1 2 5 6 9 10],:,:)
dFB = squeeze(FBgroup(:,:,1)-FBgroup(:,:,2))
avgdFB = mean(dFB,1)
dNFB = squeeze(NFBgroup(:,:,1)-NFBgroup(:,:,2))
avgdNFB = mean(dNFB,1)
% paired ttest
for T = 1:size(TF_Accuracy,2)
    [h p] = ttest(NFBgroup(:,T,1), NFBgroup(:,T,2))
    resultNFB(T,:) = [h p]
    [h p] = ttest(FBgroup(:,T,1), FBgroup(:,T,2))
    resultFB(T,:) = [h p]
    [h p] = ttest(dNFB(:,T), dFB(:,T))
    resultdFBNFB(T,:) = [h p]
end

%
[h p] = ttest2(NFBgroup(:,:,2), FBgroup(:,:,2))
%% Compare online results with calculated off line results
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];
for S = 1:12
    for TF = 1:2 % Training Feedback loop
        % online results
            if startleft(S) == 1
                ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
            else
                ear = [0 0 1 1; 1 1 0 0];
            end
         for session = 1:4
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                if length(difference) > 39
                    difference(40:end) = [];
                end
                Session_Accuracy(S,session,1) = sum(difference > 0)/length(difference)*100           
            else
                if S == 3 && session == 1 && ear(TF,session) == 2
                    load(['FB_S3_ses1_R_corrected_diff&Acc.mat'], 'difference')
                elseif S == 8 && session == 2 && ear(TF,session) == 2
                    load(['FB_S8_ses2_R_corrected_diff&Acc.mat'], 'difference')
                else
                    load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                    if length(difference) > 39
                        difference(40:end) = [];
                    end
                end
                Session_Accuracy(S,session+4,1) = sum(difference > 0)/length(difference)*100
            end
         end
    end
     % offline analysis
     cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr2')
     load(['Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr2.mat'], 'session_Accuracy')
     Session_Accuracy(S, 1:4, 2) = session_Accuracy(1, 1:4, 1)
     Session_Accuracy(S, 5:8, 2) = session_Accuracy(1, 1:4, 2)
     verschil(S,:) = Session_Accuracy(S,:,1)-Session_Accuracy(S,:,2)
end

SubA = squeeze(mean(Session_Accuracy,2));

% Significant difference? -> Yes!?
% difference is very small...
start = 1
for S = 1:12
    for ses = 1:8
        stp = start + 7;
        lange(start : stp, :) = Session_Accuracy(S,:,:)
        start = stp;
    end
end

[h p]=ttest(lange(:,1), lange(:,2))
[h2 p2] = ttest(SubA(:,1), SubA(:,2))

%% Online subject averages maken en samenvoegen
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

for S = 1:12
    for TF = 1:2 % Training Feedback loop

            if startleft(S)
                ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
            else
                ear = [0 0 1 1; 1 1 0 0];
            end

         for session = 1:4
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
            if ear(TF, session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
            if TF == 1
                load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
            else
                if S == 3 && session == 1 && ear(TF,session) == 2
                    load(['FB_S3_ses1_R_corrected_diff&Acc.mat'], 'difference')
                elseif S == 8 && session == 2 && ear(TF,session) == 2
                    load(['FB_S8_ses2_R_corrected_diff&Acc.mat'], 'difference')
                else
                    load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference') 
                end
            end
            if length(difference)> 39
                difference(40:end) = [];
            end
            Session_Accuracy(S,(session+(TF-1)*4)) = sum(difference > 0)/length(difference)*100    
         end
    end
end
Subject_Accuracy_online = mean(Session_Accuracy,2)

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save(['Online experimentes accuracies put together'], 'Session_Accuracy', 'Subject_Accuracy_online')
%% Correlation with number of correct questions? 
load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Subject_Accuracy_onr1.mat'])

x = [37 32 40 30 35 32 35 33 34 36 31 29]
x2 = [37 32 30 35 32 35 33 34 36 31 29]



for trl = 1:8  
    %subplot(3,3,trl)
    figure
    if trl == 7 % Average subject accuracy over all trl
        SAM = mean(Subject_Accuracy,2); 
        y = SAM
        X = [ones(length(x),1) x'];
        b{trl} = X\y
        [c{trl} ,p{trl}] = corrcoef([x', y])
        % S3 seems like an outlier
        y2 = SAM([1 2 4 5 6 7 8 9 10 11 12],1)
        X2 = [ones(length(x2),1) x2'];
        b2{trl} = X2\y2
        [c2{trl}, p2{trl}] = corrcoef([x2', y2])
        scatter(x, y, 'b')
        hold on
        %scatter(x2,y2, 'r')
        plot(x, X*b{trl}, 'b')
        plot(x2, X2*b2{trl}, 'r')
        hold off
        title(['Average over all trial lengths'])  %Correlation between subject accuracy and the number of correct answers')
        xlabel('Number of correct answers')
        ylabel('Subject Accuracy [%]')
        %legend(['y = ' num2str(round(b(2),2)) 'x + ' num2str(round(b(1),2)) ', c = ' num2str(round(c{trl}(1,2),2)) ', p = ' num2str(round(p{trl}(1,2),4))])

    elseif trl == 8 % Online accuracy
        load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\Online experimentes accuracies put together'], 'Subject_Accuracy_online')
        y = Subject_Accuracy_online
        b{trl} = X\y
        [c{trl} ,p{trl}] = corrcoef([x', y])
        % S3 seems like an outlier
        y2 = Subject_Accuracy_online([1 2 4 5 6 7 8 9 10 11 12],1)
        X2 = [ones(length(x2),1) x2'];
        b2{trl} = X2\y2
        [c2{trl}, p2{trl}] = corrcoef([x2', y2])
        scatter(x, y, 'b')
        hold on
        %scatter(x2,y2, 'r')
        plot(x, X*b{trl}, 'b')
        plot(x2, X2*b2{trl}, 'r')
        hold off
        title(['Online accuracy'])  %Correlation between subject accuracy and the number of correct answers')
        xlabel('Number of correct answers')
        ylabel('Subject Accuracy [%]')
        %legend(['y = ' num2str(round(b(2),2)) 'x + ' num2str(round(b(1),2)) ', c = ' num2str(round(c{trl}(1,2),2)) ', p = ' num2str(round(p{trl}(1,2),4))])

    else
        y = Subject_Accuracy(:,trl)
        X = [ones(length(x),1) x'];
        b{trl} = X\y
        [c{trl} ,p{trl}] = corrcoef([x', y])
        % S3 seems like an outlier
        y2 = Subject_Accuracy([1 2 4 5 6 7 8 9 10 11 12],trl)
        X2 = [ones(length(x2),1) x2'];
        b2{trl} = X2\y2
        [c2{trl}, p2{trl}] = corrcoef([x2', y2])
        scatter(x, y, 'b')
        hold on
        %scatter(x2,y2, 'r')
        plot(x, X*b{trl}, 'b')
        plot(x2, X2*b2{trl}, 'r')
        hold off
        title(['Trial length of ' num2str(trl*10) 's'])  %Correlation between subject accuracy and the number of correct answers')
        xlabel('Number of correct answers')
        ylabel('Subject Accuracy [%]')
        %legend(['y = ' num2str(round(b(2),2)) 'x + ' num2str(round(b(1),2)) ', c = ' num2str(round(c{trl}(1,2),2)) ', p = ' num2str(round(p{trl}(1,2),4))])
    end
end
%suptitle('Correlation between accuracy and number of correct answers')
    %legend( 'data1', 'data1', 'data2', 'data2','data3', 'data3','data4', 'data4','data5', 'data5','data6', 'data6')

% % S3 seems like an outlier
% x2 = [37 32 30 35 32 35 33 34 36 31 29]
% y2 = Subject_Accuracy([1 2 4 5 6 7 8 9 10 11 12],1)
% X2 = [ones(length(x2),1) x2'];
% b2 = X2\y2
% [c2, p2] = corrcoef([x2', y2])

B = []; C = []; P = [];
for trl = 1:8
    B = [B; [b{trl}(:,1)'; b2{trl}(:,1)']]
    C = [C; c{trl}(1,2); c2{trl}(1,2)];
    P = [P; p{trl}(1,2); p2{trl}(1,2)];
end

%% Evaluate subband method
load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr6\Subject_Accuracy_onr6.mat'])
SAsubb = Subject_Accuracy;
load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Subject_Accuracy_onr1.mat'])
SAori = Subject_Accuracy(:,[1 3 6]);

figure
subplot(1,2,1)
plot([10 30 60],SAori)
subplot(1,2,2)
plot([10 30 60], SAsubb)

% Mean
mSAsubb = mean(SAsubb,1);
sSAsubb = std(SAsubb,1);
mSAori = mean(SAori,1);
sSAori = std(SAori,1);

figure
plot([10 30 60],mSAori, 'r')
hold on
plot([10 30 60],mSAsubb, 'b')
plot([10 30 60],mSAori + sSAori, '--r')
plot([10 30 60],mSAori - sSAori, '--r')
plot([10 30 60],mSAsubb + sSAsubb, '--b')
plot([10 30 60],mSAsubb - sSAsubb, '--b')
axis([10 60 75 101])
xlabel('Trial length [s]')
ylabel('Accuracy [%]')
title('Subband envelope method')
legend('Original','Subband', )

res = []
for t = 1:3
    [h p]=ttest(SAsubb(:,t), SAori(:,t))
    res(t,:) = [h p]
end

%% Effect of high pass filtering the audio
load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr3\Subject_Accuracy_onr3.mat'])
SAaudioLP = Subject_Accuracy;
load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Subject_Accuracy_onr1.mat'])
SAori = Subject_Accuracy;

figure
subplot(1,2,1)
plot(10:10:60,SAori)
subplot(1,2,2)
plot(10:10:60, SAaudioLP)

% Mean
mSAaudioLP = mean(SAaudioLP,1);
sSAaudioLP = std(SAaudioLP,1);
mSAori = mean(SAori,1);
sSAori = std(SAori,1);

figure
plot(10:10:60,mSAori, 'r')
hold on
plot(10:10:60,mSAaudioLP, 'b')
plot(10:10:60,mSAori + sSAori, '--r')
plot(10:10:60,mSAori - sSAori, '--r')
plot(10:10:60,mSAaudioLP + sSAaudioLP, '--b')
plot(10:10:60,mSAaudioLP - sSAaudioLP, '--b')
axis([10 60 75 102])
xlabel('Trial length [s]')
ylabel('Accuracy [%]')
title('High pass filtering the audio')
legend('Original','HP filtering')
hold off

for t = 1:6
    [h p]=ttest(SAaudioLP(:,t), SAori(:,t))
    res(t,:) = [h p]
end


for S = 1:12
    figure(3)
    plot(10:10:60,SAori(S,:))
    hold on
    plot(10:10:60, SAaudioLP(S,:))
    title(['Subject ' num2str(S)])
    legend('Ori', 'LP')
    hold off
    pause
end

%% Filter values: Subject Accuracies samenvoegen
subjects = 1:12
for onr = 4
for S = subjects
    tell = 0;
    NSR = 20; 
    for LP = [0.5 1 2]
        for HP = [7 8 9]
            tell = tell +1
            cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
            load(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP' num2str(HP) '_LP' num2str(LP) '_onr' num2str(onr) '.mat'], 'subject_Accuracy')
            if length(subject_Accuracy) > 3
                Subject_Accuracy(S,tell, :) = subject_Accuracy([1 3 6])
            else
                Subject_Accuracy(S,tell, :) = subject_Accuracy
            end
        end
    end
end

cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr' num2str(onr)])
save(['Subject_Accuracy_onr' num2str(onr) '.mat'], 'Subject_Accuracy')
end

%% Filter values: Analysis
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr4\Subject_Accuracy_onr4.mat')
SAf = Subject_Accuracy;

mSAf = squeeze(mean(Subject_Accuracy,1))
sSAF = std(Subject_Accuracy,1)

plot([10 30 60],mSAf)
legend('LP 0.5 HP 7', 'LP 0.5 HP 8', 'LP 0.5 HP 9', 'LP 1 HP 7', 'LP 1 HP 8', 'LP 1 HP 9', 'LP 2 HP 7', 'LP 2 HP 8', 'LP 2 HP 9')
title('Evaluating the filter values')
xlabel('Trial length [s]')
ylabel('Accuracy [%]')

% mean over LP values
m05 = mean(mSAf(1:3,:))
m1 = mean(mSAf(4:6,:))
m2 = mean(mSAf(7:9,:))

figure
plot([10 30 60],m05)
hold on 
plot([10 30 60],m1)
plot([10 30 60],m2)
legend('0.5', '1', '2')
title('Averages for each high pass filter value')
xlabel('Trial length [s]')
ylabel('Accuracy [%]')
hold off
%significance
for tell = 1:9
    for t = 1:3
        [h p]=ttest(SAf(:,tell,t), SAf(:,5,t))
        res(tell, t) = p;
    end
end
%%
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr4\Subject_Accuracy_onr4.mat')
SAf = Subject_Accuracy;

% ANOVA

for TRL = 1:3
    Y = [];
    for l = 1:3
        Y = cat(1, Y, squeeze(SAf(:, 1+3*(l-1):3*l, TRL)))
    end
    eval(['[p' num2str(TRL) ', tbl' num2str(TRL) ', stats' num2str(TRL) '] = anova2(Y,12)'])
    
    % ttest per Low
    T05 = [squeeze(Y(1:12, 1)); squeeze(Y(1:12, 2)); squeeze(Y(1:12, 3))];
    T1 = [squeeze(Y(13:24, 1)); squeeze(Y(13:24, 2)); squeeze(Y(13:24, 3))];
    T2 = [squeeze(Y(25:end, 1)); squeeze(Y(25:end, 2)); squeeze(Y(25:end, 3))];

    [h051(TRL), p051(TRL)] = ttest(T05, T1, 'Tail', 'left')
    [h12(TRL), p12(TRL)] = ttest(T1, T2, 'Tail', 'right')
    [h052(TRL), p052(TRL)] = ttest(T05, T2)
    
    M105(TRL) = mean(T1 - T05); 
    S105(TRL) = std(T1- T05);

    M12(TRL) = mean(T1 - T2); 
    S12(TRL) = std(T1- T2);
    
    M052(TRL) = mean(T05 - T2); 
    S052(TRL) = std(T05 - T2);
end
P = [p051; p12; p052]
M = [M105; M12; M052]
S = [S105; S12; S052]
%% Lags: Put all together
nsr = [10 20 30 40]
for n = 1:4
    NSR = nsr(n)
    for S = 1:12
        cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr11')
        load(['Accuracies_S ' num2str(S) '_NSR' num2str(NSR) '_HP8_LP1_onr11.mat'], 'subject_Accuracy')
        Accuracy{n}(S,:,:) = subject_Accuracy
    end
end

save('All_Subject_Accuracies', 'Accuracy')

%% Lags: analyse seperate lags
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr11')
save('All_Subject_Accuracies', 'Accuracy')

% per trial length period
for t = 1:3 
    figure
    for n = 1:4 
        subplot(2,2,n)
        x = 0:(500/(size(Accuracy{n}, 3)-1)):500;
        plot(x, squeeze(Accuracy{n}(:,t,:))')
        hold on
        plot(x, mean(squeeze(Accuracy{n}(:,t,:)))','k', 'LineWidth', 3)
        hold off
    end
end

%averaged over the 3 trial lengths
figure
nsr = [10 20 30 40]
for n = 1:4 
    subplot(2,2,n)
    x = 0:(500/(size(Accuracy{n}, 3)-1)):500;
    plot(x, squeeze(mean(Accuracy{n},2))')
    hold on
    plot(x, mean(squeeze(mean(Accuracy{n},2)))','k', 'LineWidth', 3)
    s = squeeze(std(mean(Accuracy{n},2)))
    plot(x, mean(squeeze(mean(Accuracy{n},2)))'+s,'--k', 'LineWidth', 2)
    plot(x, mean(squeeze(mean(Accuracy{n},2)))'-s,'--k', 'LineWidth', 2)

    title(['NSR = ' num2str(nsr(n))])
    xlabel('Lag [ms]')
    ylabel('Accuracy [%]')
    hold off
end
suptitle('Accuracy calculated on seperate lags')
legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'Average')

n=4
figure
x = 0:(500/(size(Accuracy{n}, 3)-1)):500;
    plot(x, squeeze(mean(Accuracy{n},2))')
    hold on
    plot(x, mean(squeeze(mean(Accuracy{n},2)))','k', 'LineWidth', 3)
    s = squeeze(std(mean(Accuracy{n},2)))
    plot(x, mean(squeeze(mean(Accuracy{n},2)))'+s,'--k', 'LineWidth', 2)
    plot(x, mean(squeeze(mean(Accuracy{n},2)))'-s,'--k', 'LineWidth', 2)

    xlabel('Lag [ms]')
    ylabel('Accuracy [%]')
    hold off
    title('Accuracy over seperate lags')
legend('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'Average')

%% Lags: put them together: onr = 10

%% NSR

%% Difference between online and mimick online
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

for S = 1:12
    % Putting the online session accuracies together. (Based on differences)
    for TF = 1:2 % Training Feedback loop
            if startleft(S) == 1
                ear = [1 1 0 0; 0 0 1 1]; % left = 1; right = 0
            else
                ear = [0 0 1 1; 1 1 0 0];
            end

         for session = 1:4
                cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
                if ear(TF, session) == 1
                    EAR = 'left'; 
                else
                    EAR = 'right';
                end
                if TF == 1                    
                    load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                else
                    load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference') 
                    if S == 3 && session == 1
                        load(['FB_S3_ses1_R_corrected_diff&Acc'], 'difference')
                    end
                end
                if length(difference) > 39
                    difference(40:end) = [];
                end
                Session_Accuracy_online(S,session,TF) = sum(difference > 0)/length(difference)*100
         end
    end
    % Mimick online
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr2')
    load(['Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr2.mat'], 'session_Accuracy')
    Mimick_Accuracy(S,:,:) = squeeze(session_Accuracy(1,:,:))
end

% difference
diffMO = Session_Accuracy_online - Mimick_Accuracy

% significance of the difference (paired ttest)
% per subject
[h1 p1] = ttest([squeeze(Session_Accuracy_online(:,:,1)), squeeze(Session_Accuracy_online(:,:,2))]' , [squeeze(Mimick_Accuracy(:,:,1)), squeeze(Mimick_Accuracy(:,:,2))]')
% overall
AO = reshape([squeeze(Session_Accuracy_online(:,:,1)), squeeze(Session_Accuracy_online(:,:,2))], [1,12*8])
AM = reshape([squeeze(Mimick_Accuracy(:,:,1)), squeeze(Mimick_Accuracy(:,:,2))], [1,12*8])
[h2 p2] = ttest(AO - AM)