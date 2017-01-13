%% Getting the CA and CUA of the feedback
% Initial parameters

clear all
subjects = 1:12
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

lag = [100 400];
LP = 1;
HP = 8;
NSR = 20;

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; 

%% Calculating the CUA and CA for the FEEDBACK sessions and putting them together     
for S = 1:12
    if startleft(S)
        ear = [0 0 1 1]; % left = 1; right = 0
    else
        ear = [1 1 0 0];
    end
    
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])
    load(['Best_decoder_Subject' num2str(S) '.mat'], 'avg_dec')
        
    for session = 1:4 % 4 right ear session - Train on each one separate
        if S == 3 && session == 1
            load(['FB_S3_ses1_R_corrected_diff&Acc.mat'], 'CA', 'CUA', 'difference')
            difference_f(S, session, :) = difference; 
            CA_f(S, session, :) = CA; 
            CUA_f(S, session, :) = CUA;
            clear difference CA CUA
        elseif S == 8 && session == 2
            load(['FB_S8_ses2_R_corrected_diff&Acc.mat'], 'CA', 'CUA', 'difference')
            difference_f(S, session, :) = difference; 
            CA_f(S, session, :) = CA; 
            CUA_f(S, session, :) = CUA;
            clear difference CA CUA
        else
            if ear(session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end

            load(['Feedback_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'EE', 'Attended', 'Unattended', 'difference')

            if length(difference) > 39
                difference(40:end) = [];
                Attended(40:end) = []; 
                Unattended(40:end) = [];
            end
            difference_check = difference; 
            clear difference

            noftrials = size(Attended,2);
            train = (1:noftrials);

            for HH = 1:noftrials % repeat for all trials - leave one out
                trials2 = train(1:end-1);

                recenv = EE{train(end)} * avg_dec;

                CA_f(S, session, train(end)) = corr2(recenv,Attended{train(end)}); % Attended
                CUA_f(S, session, train(end)) = corr2(recenv,Unattended{train(end)}); % Unattended
                train = [train(end) train(1:end-1)];
            end       
            % check difference
            if squeeze(difference_check) == squeeze(CA_f(S, session, :) - CUA_f(S, session, :))'; 
            else
                warning('difference_checkk is not equal to new difference')
                w(S,session) = 1
            end
            difference_f(S, session, :) = squeeze(CA_f(S, session, :) - CUA_f(S, session, :)); 
        end
    end
end

difference = difference_f
CA = CA_f
CUA = CUA_f
%% save

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments')
save(['CA and CUA of online feedback.mat'], 'difference', 'CUA', 'CA')

%% Calculating the CUA and CA for the TRAINING sessions and putting them together     
clear all

subjects = 1:12
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

lag = [100 400];
LP = 1;
HP = 8;
NSR = 20;

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; 

%
for S = 1:12
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject' num2str(S)])

    if startleft(S)
        ear = [1 1 0 0]; % left = 1; right = 0
    else
        ear = [0 0 1 1];
    end
    
    
    if S < 4 % Still need to calculate CA and CUA
        for session = 1:4
            if ear(session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
    
            load(['Decoder_Params_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'])
            noftrials = 39;
            train = (1:noftrials);
            for HH = 1:noftrials % repeat for all trials - leave one out
                trials2 = train(1:end-1);

                RXX = squeeze(mean(Rxx(trials2,:,:),1)); % plain averaging of cov matrices
                RXY_ATT = squeeze(mean(Rxy_at(trials2,:),1))';

                avg_dec = RXX \ RXY_ATT; % LS solution.
                recenv = EE{train(end)} * avg_dec;

                CA_t(S, session, train(end)) = corr2(recenv,Attended{train(end)}); % Attended
                CUA_t(S, session, train(end)) = corr2(recenv,Unattended{train(end)}); % Unattended
                train = [train(end) train(1:end-1)];
            end

            difference_t(S, session, :) = CA_t(S, session, :) - CUA_t(S, session, :);

        end
    else % CUA and CA were calculated here...
        for session = 1:4
            if ear(session) == 1
                EAR = 'left'; 
            else
                EAR = 'right';
            end
    
            load(['Session_Accuracy_Subject' num2str(S) '_Session' num2str(session) '_' EAR '.mat'], 'CA', 'CUA', 'difference')
            CA_t(S, session, :) = CA;
            CUA_t(S, session, :) = CUA;
            difference_t(S, session, :) = difference;
            
            clear CA CUA difference
        end
            
    end
end
        
%% save
difference = difference_t
CA = CA_t
CUA = CUA_t

cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments')
save(['CA and CUA of online training.mat'], 'difference', 'CUA', 'CA')


%% Split in FB and NFB groups (NIET GEBRUIKT...)
% ONLINE
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments')
load(['CA and CUA of online training.mat'], 'difference', 'CUA', 'CA')
difference_FB(1,:,:,:) = difference([1 2 5 6 9 10], :, :);
difference_NFB(1,:,:,:) = difference([3 4 7 8 11 12], :, :);
CA_FB(1,:,:,:) = CA([1 2 5 6 9 10], :, :);
CA_NFB(1,:,:,:) = CA([3 4 7 8 11 12], :, :);
CAU_FB(1,:,:,:) = CAU([1 2 5 6 9 10], :, :);
CAU_NFB(1,:,:,:) = CAU([3 4 7 8 11 12], :, :);

load(['CA and CUA of online feedback.mat'], 'difference', 'CUA', 'CA')
difference_FB(1,:,:,:) = difference([1 2 5 6 9 10], :, :);
difference_NFB(1,:,:,:) = difference([3 4 7 8 11 12], :, :);
CA_FB(2,:,:,:) = CA([1 2 5 6 9 10], :, :);
CA_NFB(2,:,:,:) = CA([3 4 7 8 11 12], :, :);
CAU_FB(2,:,:,:) = CAU([1 2 5 6 9 10], :, :);
CAU_NFB(2,:,:,:) = CAU([3 4 7 8 11 12], :, :);


%% Try n way ANOVA -> werkt niet... WAAROM??? 
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments')
load(['CA and CUA of online training.mat'], 'difference', 'CUA', 'CA')
CA_t = CA; CUA_t = CUA; difference_t = difference; 
clear CA CUA difference
load(['CA and CUA of online feedback.mat'], 'difference', 'CUA', 'CA')
CA_f = CA; CUA_f = CUA; difference_f = difference;
clear CA CUA difference

CA_all = [];
for S = 1:12
    for session = 1:4
        CA_all = [CA_all, squeeze(CA_t(S, session,:))', squeeze(CA_f(S,session,:))']
    end
end

g1 = [ones(1,312), ones(1,312)*2, ones(1,312)*3, ones(1,312)*4,ones(1,312)*5, ones(1,312)*6, ones(1,312)*7, ones(1,312)*8,ones(1,312)*9, ones(1,312)*10, ones(1,312)*11, ones(1,312)*12]; % Subjects
for i = 1:24
    if mod(i,2)==0
        g2(156*(i-1)+1 : i*156) = {'F'}
    else
        g2(156*(i-1)+1 : i*156) = {'T'}
    end
end

for i = 1:6
    if mod(i,2)==0
        g3(624*(i-1)+1 : i*624) = {'NFB'}
    else
        g3(624*(i-1)+1 : i*624) = {'FB'}
    end
end

[p,tbl,stats,terms] = anovan(CA_all, {g1; g2; g3}, 'model', 'interaction', 'varnames',{'g1','g2','g3'})

%% Make Onr1 to calculate beneath
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
for S = 1:12
    load(['Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'], 'subject_CA', 'subject_CUA', 'subject_difference')
    difference_t(S,:) = subject_difference{1}(1:156);
    difference_f(S,:) = subject_difference{1}(157:end);
    
    CA_t(S,:) = subject_CA{1}(1:156);
    CA_f(S,:) = subject_CA{1}(157:end);
    
    CUA_t(S,:) = subject_CUA{1}(1:156);
    CUA_f(S,:) = subject_CUA{1}(157:end);
end

%% Verschil tussen training en feedback sessies per subject: t-test
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments')
load(['CA and CUA of online training.mat'], 'difference', 'CUA', 'CA')

difference_t = []; CUA_t = []; CA_t = [];

for session = 1:4
    difference_t = [difference_t, squeeze(difference(:,session,:))]; 
    CUA_t = [CUA_t, squeeze(CUA(:,session,:))]; 
    CA_t = [CA_t, squeeze(CA(:,session,:))];;
end

clear difference CA CUA

load(['CA and CUA of online feedback.mat'], 'difference', 'CUA', 'CA')
difference_f = []; CUA_f = []; CA_f = [];

for session = 1:4
    difference_f = [difference_f, squeeze(difference(:,session,:))]; 
    CUA_f = [CUA_f, squeeze(CUA(:,session,:))]; 
    CA_f = [CA_f, squeeze(CA(:,session,:))];;
end

clear difference CA CUA
%%

for S = 1:12
    % test normality
    [hndiff_t(S, 1), pndiff_t(S, 1)] = adtest(difference_t(S,:))
    [hndiff_f(S, 1), pndiff_f(S, 1)] = adtest(difference_f(S,:))
    [hndiff_t(S, 2), pndiff_t(S, 2)] = chi2gof(difference_t(S,:))
    [hndiff_f(S, 2), pndiff_f(S, 2)] = chi2gof(difference_f(S,:))
    [hndiff_t(S, 3), pndiff_t(S, 3)] = jbtest(difference_t(S,:))
    [hndiff_f(S, 3), pndiff_f(S, 3)] = jbtest(difference_f(S,:))
    
    [hnCA_t(S, 1), pnCA_t(S, 1)] = adtest(difference_t(S,:))
    [hnCA_f(S, 1), pnCA_f(S, 1)] = adtest(difference_f(S,:))
    [hnCA_t(S, 2), pnCA_t(S, 2)] = chi2gof(difference_t(S,:))
    [hnCA_f(S, 2), pnCA_f(S, 2)] = chi2gof(difference_f(S,:))
    [hnCA_t(S, 3), pnCA_t(S, 3)] = jbtest(difference_t(S,:))
    [hnCA_f(S, 3), pnCA_f(S, 3)] = jbtest(difference_f(S,:))
    
    [hnCUA_t(S, 1), pnCUA_t(S, 1)] = adtest(difference_t(S,:))
    [hnCUA_f(S, 1), pnCUA_f(S, 1)] = adtest(difference_f(S,:))
    [hnCUA_t(S, 2), pnCUA_t(S, 2)] = chi2gof(difference_t(S,:))
    [hnCUA_f(S, 2), pnCUA_f(S, 2)] = chi2gof(difference_f(S,:))
    [hnCUA_t(S, 3), pnCUA_t(S, 3)] = jbtest(difference_t(S,:))
    [hnCUA_f(S, 3), pnCUA_f(S, 3)] = jbtest(difference_f(S,:))

    %ttest
    [hdiff(S) pdiff(S)] = ttest2(squeeze(difference_t(S,:)), squeeze(difference_f(S,:)))
    [hCUA(S) pCUA(S)] = ttest2(squeeze(CUA_t(S,:)), squeeze(CUA_f(S,:)))
    [hCA(S) pCA(S)] = ttest2(squeeze(CA_t(S,:)), squeeze(CA_f(S,:)))    
end

Mdiff = [squeeze(mean(difference_t,2))'; squeeze(mean(difference_f,2))']
Sdiff = [squeeze(std(difference_t, [],2))'; squeeze(std(difference_f,[],2))']

MCUA = [squeeze(mean(CUA_t,2))'; squeeze(mean(CUA_f,2))']
SCUA = [squeeze(std(CUA_t, [],2))'; squeeze(std(CUA_f,[],2))']

MCA = [squeeze(mean(CA_t,2))'; squeeze(mean(CA_f,2))']
SCA = [squeeze(std(CA_t, [],2))'; squeeze(std(CA_f,[],2))']

dMdiff_FB = (Mdiff(1,[1 2 5 6 9 10]) - Mdiff(2,[1 2 5 6 9 10]));
dMdiff_NFB = (Mdiff(1,[3 4 7 8 11 12]) - Mdiff(2,[3 4 7 8 11 12]));
dMCA_FB = (MCA(1,[1 2 5 6 9 10]) - MCA(2,[1 2 5 6 9 10])); 
dMCA_NFB = MCA(1,[3 4 7 8 11 12]) - MCA(2,[3 4 7 8 11 12]);
dMCUA_FB = (MCUA(1,[1 2 5 6 9 10]) - MCUA(2,[1 2 5 6 9 10])); 
dMCUA_NFB = MCUA(1,[3 4 7 8 11 12]) - MCUA(2,[3 4 7 8 11 12]);

% normality test OPM: not enough data for CHI squared
    [hndiff_FB(1), pndiff_FB(1)] = adtest(dMdiff_FB)
    [hndiff_NFB(1), pndiff_NFB(1)] = adtest(dMdiff_NFB)
%     [hndiff_FB(2), pndiff_FB(2)] = chi2gof(dMdiff_FB)
%     [hndiff_NFB(2), pndiff_NFB(2)] = chi2gof(dMdiff_NFB)
    [hndiff_FB(3), pndiff_FB(3)] = jbtest(dMdiff_FB)
    [hndiff_NFB(3), pndiff_NFB(3)] = jbtest(dMdiff_NFB)
    
    [hnCA_FB(1), pnCA_FB(1)] = adtest(dMCA_FB)
    [hnCA_NFB(1), pnCA_NFB(1)] = adtest(dMCA_NFB)
%     [hnCA_FB(2), pnCA_FB(2)] = chi2gof(dMCA_FB)
%     [hnCA_NFB(2), pnCA_NFB(2)] = chi2gof(dMCA_NFB)
    [hnCA_FB(3), pnCA_FB(3)] = jbtest(dMCA_FB)
    [hnCA_NFB(3), pnCA_NFB(3)] = jbtest(dMCA_NFB)
    
    [hnCUA_FB(1), pnCUA_FB(1)] = adtest(dMCUA_FB)
    [hnCUA_NFB(1), pnCUA_NFB(1)] = adtest(dMCUA_NFB)
%     [hnCUA_FB(2), pnCUA_FB(2)] = chi2gof(dMCUA_FB)
%     [hnCUA_NFB(2), pnCUA_NFB(2)] = chi2gof(dMCUA_NFB)
    [hnCUA_FB(3), pnCUA_FB(3)] = jbtest(dMCUA_FB)
    [hnCUA_NFB(3), pnCUA_NFB(3)] = jbtest(dMCUA_NFB)

% ttest
[hMdiff pMdiff] = ttest2(dMdiff_FB, dMdiff_NFB)
[hMCA pMCA] = ttest2(dMCA_FB, dMCA_NFB)
[hMCA pMCUA] = ttest2(dMCUA_FB, dMCUA_NFB)
    
