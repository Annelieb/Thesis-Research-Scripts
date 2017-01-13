%% MAKE & TEST DECODERS

close all 
clear all 

%% Initial parameters
nofsubj = 4;

%% MAKE SESSION DECODERS

for TRL = 1:6
    cd('C:\Users\Annelies\Documents\Thesis\Blijft op one drive\Thesis one drive\script\All_data\New bp filter\intermediate\Intermediate 500')
    load(['TRL' num2str(TRL) '_ISR8_NSR3.mat'])
    
    for S = 1:nofsubj
        sessionlength = size(covar.Rxy_at_all,2)/8;
        session = [];
        for ses = 1:8
            avg_dec = []; RXX = []; RXY_ATT = [];
            session = ((1+ sessionlength *(ses-1)):(sessionlength*ses));
            RXX = squeeze(mean(covar.Rxx_all(S,session,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,session,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            
            cd('C:\Users\Annelies\Documents\Thesis\Blijft op one drive\Thesis one drive\script\All_data\Decoders')
            save(['Decoder_Subj' num2str(S) '_Ses' num2str(ses) '_TRL' num2str(TRL) '.mat'], 'avg_dec')
        end
    end
end

%% TEST SESSION DECODERS

for TRL = 1:6
    cd('C:\Users\Annelies\Documents\Thesis\Blijft op one drive\Thesis one drive\script\All_data\New bp filter\intermediate\Intermediate 500')
    load(['TRL' num2str(TRL) '_ISR8_NSR3.mat'])
    
    for S = 1:nofsubj
        noftrials = size(covar.Rxy_at_all,2);
        sessionlength = noftrials/8;
        nofothers = noftrials-sessionlength;
        trials = (1:noftrials);
        

        for dec = 1:8
            cd('C:\Users\Annelies\Documents\Thesis\Blijft op one drive\Thesis one drive\script\All_data\Decoders')
            load(['Decoder_Subj' num2str(S) '_Ses' num2str(dec) '_TRL' num2str(TRL) '.mat'])
            
            others = trials(1+sessionlength:end);
            CA = zeros(1,noftrials); CUA = zeros(1,noftrials);
            
            for HH = 1:nofothers                
                recenv = EE{S}(:,:,others(HH))*avg_dec;
                
                CA(others(HH)) = corr2(recenv,Attended{S}(others(HH),:)');
                CUA(others(HH)) = corr2(recenv,Unattended{S}(others(HH),:)');
            end
            trials = [trials(1+sessionlength:end),trials(1:sessionlength)];
            
            Accuracy(TRL, S, dec) =  (sum(CA-CUA> 0)/(nofothers))*100
            
            if dec < 5
                Accuracy_left(TRL, S, dec) = (sum((CA(1:noftrials/2)-CUA(1:noftrials/2))> 0)/(3*sessionlength))*100
                Accuracy_right(TRL, S, dec) = (sum((CA(noftrials/2+1 : end)-CUA(noftrials/2+1 : end))> 0)/(4*sessionlength))*100
            elseif dec > 4
                Accuracy_left(TRL, S, dec) = (sum((CA(1:noftrials/2)-CUA(1:noftrials/2))> 0)/(4*sessionlength))*100
                Accuracy_right(TRL, S, dec) = (sum((CA(noftrials/2+1 : end)-CUA(noftrials/2+1 : end))> 0)/(3*sessionlength))*100
            end

        end
    end
end

cd('C:\Users\Annelies\Documents\Thesis\Blijft op one drive\Thesis one drive\script\All_data\Decoders')
save(['Result_session_decoders.mat'], 'Accuracy', 'Accuracy_left', 'Accuracy_right')



