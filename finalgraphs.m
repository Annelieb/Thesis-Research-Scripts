load('C:\Users\Annelies\Documents\Studie\Thesis\script\Before and after realizing mistake in script\Results_A.mat')
correct = Subject_Accuracy;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\Before and after realizing mistake in script\Results_R.mat')
wrong = Subject_Accuracy;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\ori\New_bp_ori.mat')
newbp = Accuracy_new_bp;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate500.mat')
int = Intermediate500;

mc = mean(correct,2);
mw = mean(wrong,2);
mnbp = mean(newbp,2);
mint = mean(int,2);

% Significance level
noftrials = [312, 152, 104, 72, 56, 48];

for TRL = 1:6
    N = noftrials(TRL); 
    siglevel(TRL) = binoinv(0.95, N, 0.5)/N*100;
end


% Comparing the averaged accuracies
figure
plot(10:10:60,mw)
hold on
plot(10:10:60, mc)
plot(10:10:60,mnbp)
plot(10:10:60,mint)
plot(10:10:60, siglevel, ':k')

legend('1. Wrong left out', '2. Correct left out', '3. New bp filter', '4. Down sampling after enveloping', 'Chance level', 'Location', 'southeast')
title('A step by step comparison of the algorithm improvement')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

not2 = mean(int(:,[1 3 4]),2);

% Comparing the final accuracies
figure
plot(10:10:60, int)
hold on
plot(10:10:60, mint, '--k')
%plot(10:10:60, not2, '-.k')
plot(10:10:60, siglevel, ':k')
legend('Subject 1', 'Subject 2', 'Subject 3', 'Subject 4', 'Average', 'Chance level', 'Location', 'southeast')
title('Subject accuracies after the error correction')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% General decoders
% Make the general decoders (Because I cannot find where I had made these
% previously)
for TRL = 1:6
    load(['C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate 500\TRL' num2str(TRL) '_ISR8_NSR3.mat'])
    trainS = [1,3,4]; % Decide here whether or not subject 2 is included.
    Rxx_all = []; Rxy_all = [];
    
    for S = trainS 
        trialsS = trainS(2:end)
        
        for tr = 1:length(trialsS)
            Rxx_all = cat(1,Rxx_all,squeeze(covar.Rxx_all(trialsS(tr),:,:,:)));
            Rxy_all = cat(1,Rxy_all,squeeze(covar.Rxy_at_all(trialsS(tr),:,:)));
        end
        
        RXX = squeeze(mean(Rxx_all,1));
        RXY_ATT = squeeze(mean(Rxy_all,1));
        avg_dec = RXX \ RXY_ATT';
        
        trainS = [trainS(2:end), trainS(1)];
        
        noftrials = size(Attended{S},1);
          
        for HH = 1:noftrials % repeat for all trials - leave one out
            
            recenv = EE{S}(:,:,HH) * avg_dec;
        
            CA(HH) = corr2(recenv,Attended{S}(HH,:)'); % Attended
            CUA(HH) = corr2(recenv,Unattended{S}(HH,:)'); % Unattended
            
        end
        difference(S,:) = CA - CUA;
        
        Subject_Accuracy(TRL,S) = (sum(CA-CUA> 0)/(noftrials))*100; % Final accuracy on all trials
        CA = []; CUA = []; difference = [];
        
    end
end

% save
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\Cross subjects')
save('Cross_Subjects_not2', 'Subject_Accuracy')

%% Graphs of general decoders
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\Cross subjects')
load('Cross_Subjects_not2', 'Subject_Accuracy')
Not2 = Subject_Accuracy;
load('Cross_Subjects_all', 'Subject_Accuracy')
all = Subject_Accuracy;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate500.mat')
int = Intermediate500;


m2 = mean(Not2,2);
mall = mean(all,2);
mint = mean(int,2);

figure
plot(10:10:60, mall)
hold on
plot(10:10:60, m2)
plot(10:10:60, mint, 'k')
plot(10:10:60, siglevel, ':k')
legend('Cross-Subject decoder', 'Not Subject 2', 'Subject Specific', 'Chance level', 'Location', 'southeast')
title('Cross-Subject decoders ')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% Subband
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate500.mat')
int = Intermediate500;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\subband\Accuracy.mat')
sub = Accuracy;

mint = mean(int,2);
msub = mean(sub,2);

figure
plot(10:10:60, mint)
hold on 
plot(10:10:60, msub)
legend('No sub-band', 'Sub-band', 'Location', 'southeast')
title('Sub-band method')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% High pass filtering the audio
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate500.mat')
int = Intermediate500;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Result_LPaudiois0.mat')
Hp0 = squeeze(Subject_Accuracy(:,end,end,:));

mint = mean(int,2);
mHp0 = mean(Hp0,2);

figure
plot(10:10:60, mint)
hold on 
plot(10:10:60, mHp0)
legend('Audio HP-filtered at 2Hz', 'Audio not HP-filtered', 'Location', 'southeast')
title('High Pass (HP) filtering the audio envelope')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% Filter values
LP = [0.5, 1, 2];
HP = [7 8 9];
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Investigate LP and HP values')
for l = 1:3
    for h = 1:3
        load(['Result_LP' num2str(LP(l)) '_HP' num2str(HP(h)) '.mat'])
        SA(:,:,l,h) = squeeze(Subject_Accuracy(:,8,3,:));
    end
end

% graphs for the averages over the h (so for each l)
mSA = squeeze(mean(SA, 2));
mmSA = squeeze(mean(mSA,3));

figure
plot(10:10:60, mmSA)
legend('0.5 Hz', '1 Hz', '2 Hz', 'Location', 'southeast')
title('Comparing filter values; focus on the lower frequency limit')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% Final comparisson of the pilot study

load('C:\Users\Annelies\Documents\Studie\Thesis\script\Before and after realizing mistake in script\Results_A.mat')
correct = Subject_Accuracy;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\ori\New_bp_ori.mat')
newbp = Accuracy_new_bp;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate500.mat')
int = Intermediate500;
load('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Result_LPaudiois0.mat')
Hp0 = squeeze(Subject_Accuracy(:,end,end,:));
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Investigate LP and HP values')
load(['Result_LP' num2str(1) '_HP' num2str(8) '.mat'])
filt = squeeze(Subject_Accuracy(:,8,3,:));

mHp0 = mean(Hp0,2);
mc = mean(correct,2);
mnbp = mean(newbp,2);
mint = mean(int,2);
mfilt = mean(filt,2);
mFn2 = mean(filt(:,[1 3 4]),2);

noftrials = [312, 152, 104, 72, 56, 48];

for TRL = 1:6
    N = noftrials(TRL); 
    siglevel(TRL) = binoinv(0.95, N, 0.5)/N*100;
end

% Comparing the averaged accuracies
figure
plot(10:10:60, mc)
hold on
%plot(10:10:60,mnbp)
%plot(10:10:60,mint)
%plot(10:10:60, mHp0)
plot(10:10:60, mfilt)
plot(10:10:60, mFn2)
plot(10:10:60, siglevel, ':k')

%legend('2. Correct left out', '3. New bp filter', '4. Down sampling after enveloping', '5. Only low pass filtering the audio', '6. Filter values 1-8Hz', '7. Average without S2', 'Significance level', 'Location', 'southeast')
legend('Before improvements', 'After improvements', 'After improvements without S2', 'Significance level', 'Location', 'southeast')
title('Improvements during pilot study')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% Comparison 'cut first', 'start online', 'cut first & start online' and
% 'mimick online'
for S = 1:12
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr7\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr7.mat'], 'subject_Accuracy')
    SA_cut_first(S,:) = subject_Accuracy;
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr8\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr8.mat'], 'subject_Accuracy')
    SA_start_online(S,:) = subject_Accuracy;
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr9\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr9.mat'], 'subject_Accuracy')
    SA_cf_so(S,:) = subject_Accuracy;
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'], 'subject_Accuracy')
    SA_norm(S,:) = subject_Accuracy;
   % load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr2\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr2.mat'], 'session_Accuracy')
   % SA_mimick(S,:) = squeeze(mean(squeeze(mean(session_Accuracy,2)),2));
end

%% Summary of on-line results
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\Online experimentes accuracies put together.mat')

boxplot(Session_Accuracy')
title('Box plots of the on-line results')
xlabel('Subject')
ylabel('Decoding Accuracy [%]')

%% Same for mimmick and onr1
for S = 1:12
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr2\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr2.mat'], 'session_Accuracy')
    SA_mimick(S,:) = [squeeze(session_Accuracy(1,:,1)), squeeze(session_Accuracy(1,:,2))];
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'], 'session_Accuracy')
    SA_onr1(S,:) = [squeeze(session_Accuracy(1,:,1)), squeeze(session_Accuracy(1,:,2))];  
end

load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\Online experimentes accuracies put together.mat')

boxplot(Session_Accuracy')
hold on
boxplot(SA_mimick')
boxplot(SA_onr1')
title('Box plots of the on-line results')
legend('Online', 'Mimick', 'Offline')
xlabel('Subject')
ylabel('Decoding Accuracy [%]')

%% ANOVA of the on-line experiments... 
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\Online experimentes accuracies put together.mat')

Feedback = Session_Accuracy([1 2 5 6 9 10],:);
NoFB = Session_Accuracy([3 4 7 8 11 12],:);

mFBT = squeeze(mean(Feedback(:,1:4),2));
mFBO = squeeze(mean(Feedback(:,5:8),2));

mNoFBT = squeeze(mean(NoFB(:,1:4),2));
mNoFBO = squeeze(mean(NoFB(:,5:8),2));

% ANOVA
Y = [mNoFBT, mFBT; mNoFBO, mFBO]

[p,tbl] = anova2(Y,6)

%% ANOVA of the off-line experiments... 
for S = 1:12
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'], 'subject_difference')
    for TRL = 1:6
        noftrl = length(subject_difference{TRL});
        First4{TRL}(S) = sum(subject_difference{TRL}(1:noftrl/2) > 0)/(noftrl/2); 
        Last4{TRL}(S) = sum(subject_difference{TRL}(noftrl/2+1:end) > 0)/(noftrl/2); 
    end
end

for TRL = 1:6
    Y{TRL} = [First4{TRL}([3, 4, 7, 8, 11, 12])', First4{TRL}([1, 2, 5, 6, 9, 10])';  Last4{TRL}([3, 4, 7, 8, 11, 12])', Last4{TRL}([1, 2, 5, 6, 9, 10])'];
    figure
    eval(['[p_' num2str(TRL) ',tbl_' num2str(TRL) '] = anova2(Y{' num2str(TRL) '},6)'])
end

%% ANOVA of difference... 
% 2 factor ANOVA per subject
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];
T = []; F = [];
for S = 1:12
    T = []; F = [];
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    for TF = 1:2
        if startleft(S) == 1
            ear = [1 1 0 0; 0 0 1 1];
        else
            ear = [0 0 1 1; 1 1 0 0];
        end
        for session = 1:4
            if ear(TF, session) == 1
                EAR = 'left';
            else
                EAR = 'right';
            end
            if TF == 1
                load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                if length(difference) > 39
                    difference(40:end) = [];
                end
                T = [T, difference];
                clear difference
            else
                if S == 3 && session == 1
                    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\FB_S3_ses1_R_corrected_diff&Acc.mat'], 'difference')
                elseif S == 8 && session == 2
                    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\FB_S8_ses2_R_corrected_diff&Acc.mat'], 'difference')
                else
                    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S) '\Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'difference')
                end
                if length(difference) > 39
                    difference(40:end) = [];
                end
                F = [F, difference];
                clear difference
            end
        end
    end
    Training(S,:) = T;
    Feedback(S,:) = F;
end

% ttest for each subject
for S = 1:12
    [h(S), p(S)] = ttest2(Training(S,:), Feedback(S,:));
end

% 3 way ANOVA -> werkt niet om een of andere reden...
Y = [];
for S = 1:12
    for TF = 1:2
        if TF == 1
            Y = [Y, Training(S,:)];
        else
            Y = [Y, Feedback(S,:)];
        end
    end
end
% subjects
g1 = [ones(1,312), ones(1,312)*2, ones(1,312)*3, ones(1,312)*4, ones(1,312)*5, ones(1,312)*6, ones(1,312)*7, ones(1,312)*8, ones(1,312)*9, ones(1,312)*10, ones(1,312)*11, ones(1,312)*12]; % Subject
% TO
g2 = [ones(1,156), ones(1,156)*2, ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2, ones(1,156), ones(1,156)*2, ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2,ones(1,156), ones(1,156)*2];
% FB NFB
g3 = [ones(1,312), ones(1,312), ones(1,312)*2, ones(1,312)*2, ones(1,312), ones(1,312), ones(1,312)*2, ones(1,312)*2, ones(1,312), ones(1,312), ones(1,312)*2, ones(1,312)*2]; % Subject

p = anovan(Y, {g1, g2, g3}, 'model', 'interaction', 'varnames', {'g1', 'g2', 'g3'})

% ANOVA of averages
MF = mean(Feedback, 2)
MT = mean(Training, 2)



    Y = [MT([3, 4, 7, 8, 11, 12])', MT([1, 2, 5, 6, 9, 10])';  MF([3, 4, 7, 8, 11, 12])', MF([1, 2, 5, 6, 9, 10])'];
    figure
    [p, tbl] = anova2(Y,6)
%% For the offline onr1
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
for S = 1:12
    load(['Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'], 'subject_difference')
    for TRL = 1:6
        Training{TRL}(S,:) = subject_difference{TRL}(1:end/2);
        Online{TRL}(S,:) = subject_difference{TRL}(end/2+1:end);
    end
end

for TRL = 1:6
    for S = 1:12
        [h(S,TRL), p(S,TRL)] = ttest2(Training{TRL}(S,:), Online{TRL}(S,:))
    end
end

for TRL = 1:6
    MF = mean(Online{TRL}, 2)
    MT = mean(Training{TRL}, 2)
    Y{TRL} = [MT([3, 4, 7, 8, 11, 12]), MT([1, 2, 5, 6, 9, 10]);  MF([3, 4, 7, 8, 11, 12]), MF([1, 2, 5, 6, 9, 10])];
    figure
    eval(['[p_' num2str(TRL) ',tbl_' num2str(TRL) '] = anova2(Y{' num2str(TRL) '},6)'])
end


%% Channel reduction

for S = 1:12
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis'])
    load(['Channel_reduction2_S ' num2str(S) '_trl10.mat'], 'subject_Accuracy', 'removedchannels')
    SA(S,:) = subject_Accuracy;
    RC(S,:) = removedchannels{1,1};
end
 C = ones(1,24)*24 - [0:23];
 
MSA = mean(SA,1);
sSA = std(SA, [], 1);
figure 
plot(C, SA)
hold on
plot(C, ones(1,24)*binoinv(0.95, 312, 0.5)/312*100, ':k')
legend('S1', 'S2','S3', 'S4','S5', 'S6','S7', 'S8','S9', 'S10','S11', 'S12', 'Significance level', 'Location', 'southwest')
title('Elektrode reduction')
xlabel('Number of elektrodes')
ylabel('Accuracy [%]')
axis([1 24 45 90])
set ( gca, 'xdir', 'reverse' )

figure
plot(C, MSA, 'k')
hold on
plot(C, MSA + sSA, '--k')
plot(C, MSA - sSA, '--k')
plot(C, ones(1,24)*binoinv(0.95, 312, 0.5)/312*100, ':k')
hold off
title('Average elektrode reduction')
xlabel('Number of elektrodes')
ylabel('Decoding Accuracy [%]')
axis([1 24 45 90])
set ( gca, 'xdir', 'reverse' )

cmap = hsv(13)

figure
for S = 1:12
    plot(C, SA(S,:), 'Color', cmap(S,:))
    hold on
end
plot(C, MSA, 'k', 'LineWidth', 2)
plot(C, ones(1,24)*binoinv(0.95, 312, 0.5)/312*100, ':k')
plot(C, MSA + sSA, '--k', 'LineWidth', 2)
plot(C, MSA - sSA, '--k', 'LineWidth', 2)
hold off
legend('S1', 'S2','S3', 'S4','S5', 'S6','S7', 'S8','S9', 'S10','S11', 'S12', 'Average', 'Significance level', 'Location', 'southwest')
title('Elektrode reduction')
xlabel('Number of elektrodes')
ylabel('Decoding Accuracy [%]')
axis([1 24 45 90])
set ( gca, 'xdir', 'reverse' )

%% CR all same
load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\CR_allS_trl_1.mat'], 'subject_Accuracy')
SA = subject_Accuracy{1,1};
MSA = mean(subject_Accuracy{1,1});
sSA = std(subject_Accuracy{1,1});
cmap = hsv(13)

 C = ones(1,24)*24 - [0:23];
figure
for S = 1:12
    plot(C, SA(S,:), 'Color', cmap(S,:))
    hold on
end
cmap = hsv(13)
plot(C, MSA, 'k', 'LineWidth', 2)
plot(C, ones(1,24)*binoinv(0.95, 312, 0.5)/312*100, ':k')
plot(C, MSA + sSA, '--k', 'LineWidth', 2)
plot(C, MSA - sSA, '--k', 'LineWidth', 2)
hold off
legend('S1', 'S2','S3', 'S4','S5', 'S6','S7', 'S8','S9', 'S10','S11', 'S12', 'Average', 'Significance level', 'Location', 'southwest')
title('Elektrode reduction')
xlabel('Number of elektrodes')
ylabel('Decoding Accuracy [%]')
axis([1 24 45 90])
set ( gca, 'xdir', 'reverse' )


%% Compare Final results with Pilot results
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Investigate LP and HP values')
load(['Result_LP' num2str(1) '_HP' num2str(8) '.mat'])
pilot = squeeze(Subject_Accuracy(:,8,3,:));
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Subject_Accuracy_onr1.mat')
final = Subject_Accuracy;

mP = mean(pilot,2);
mPn2 = mean(pilot(:,[1 3 4]),2);
mF = mean(final,1);
sF = std(final);

figure
plot(10:10:60, mP)
hold on
plot(10:10:60, mPn2)
plot(10:10:60, mF, 'k')
plot(10:10:60, mF + sF, '--k')
plot(10:10:60, mF - sF, '--k')
legend('Pilot', 'Pilot without 2', 'Final', 'Location', 'southeast')
title('Comparing the pilot study and final experiment')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

%% Decoder choice
cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1'])
load(['FB_compare_decoders.mat'], 'FB_Accuracy')

% Significance level
noftrials = [312, 152, 104, 72, 56, 48]/2;

for TRL = 1:6
    N = noftrials(TRL); 
    siglevel(TRL) = binoinv(0.95, N, 0.5)/N*100;
end

%

M = squeeze(mean(FB_Accuracy,1));
S = squeeze(std(FB_Accuracy,[],1));

figure
plot(10:10:60,M)
hold on
plot(10:10:60, siglevel, ':k')
legend('1. Chosen decoder', '2. Leave one out', '3. All 4 training', '4. All 4 pilot', '5. Pilot 1, 3 and 4', '6. All other final', 'Chance level', 'Location', 'eastoutside')
title('Comparing decoder options')
ax = gca;
ax.XTick = 10:10:60
xlabel('Trial Length [s]')
ylabel('Decoding Accuracy [%]')
hold off

% Significance
for TRL = 1:6
    [h12(TRL), p12(TRL)] = ttest(squeeze(FB_Accuracy(:,TRL,1)),squeeze(FB_Accuracy(:,TRL,2)), 'Tail', 'left');
    [h13(TRL), p13(TRL)] = ttest(squeeze(FB_Accuracy(:,TRL,1)),squeeze(FB_Accuracy(:,TRL,3)), 'Tail', 'left');
    [h23(TRL), p23(TRL)] = ttest(squeeze(FB_Accuracy(:,TRL,2)),squeeze(FB_Accuracy(:,TRL,3)), 'Tail', 'right')
end
