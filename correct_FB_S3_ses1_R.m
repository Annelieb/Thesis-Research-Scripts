%% FB S3 ses1 R
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject3\Feedback_Subject3_Session1_right.mat', 'EE', 'Attended', 'Unattended')
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject3\Best_decoder_Subject3.mat', 'avg_dec')

noftrials = size(Attended,2);
train = (1:noftrials);


for HH = 1:noftrials % repeat for all trials - leave one out
    trials2 = train(1:end-1);
        
    recenv = EE{train(end)} * avg_dec;
        
    CA(train(end)) = corr2(recenv,Attended{train(end)}); % Attended
    CUA(train(end)) = corr2(recenv,Unattended{train(end)}); % Unattended
    train = [train(end) train(1:end-1)];
end

difference = CA - CUA;
Session_Accuracy = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials

save(['FB_S3_ses1_R_corrected_diff&Acc'], 'difference', 'Session_Accuracy', 'CA', 'CUA')

%% FB S8 ses2

load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject8\Feedback_Subject8_Session2_right.mat', 'EE', 'Attended', 'Unattended')
load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Subject8\Best_decoder_Subject8.mat', 'avg_dec')

noftrials = size(Attended,2);
train = (1:noftrials);

for HH = 1:noftrials % repeat for all trials - leave one out
    trials2 = train(1:end-1);
        
    recenv = EE{train(end)} * avg_dec;
        
    CA(train(end)) = corr2(recenv,Attended{train(end)}); % Attended
    CUA(train(end)) = corr2(recenv,Unattended{train(end)}); % Unattended
    train = [train(end) train(1:end-1)];
end

difference = CA - CUA;
Session_Accuracy = (sum(CA-CUA> 0)/(noftrials))*100 % Final accuracy on all trials

save(['FB_S8_ses2_R_corrected_diff&Acc'], 'difference', 'Session_Accuracy', 'CA', 'CUA')
