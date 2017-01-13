%% MAKE THE DECODERS
clear all 

% Adjustable parameters
%nofsessions =4; % Number of sessions you want the decoder to be built on. 
nofsubj = 4;

% Making the decoders
for nofsessions = 1:4;
combinations = nchoosek(1:8, nofsessions);
nofcomb = size(combinations, 1); % Total number of combinations
    
for TRL = 7:8
    cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate 500')
    load(['TRL' num2str(TRL) '_ISR8_NSR3.mat'])
    

    
    for S = 1:nofsubj
        sessionlength = size(covar.Rxy_at_all,2)/8;       
        for comb = 1:nofcomb
            avg_dec = []; RXX = []; RXY_ATT = [];
            sessions = [];
            
            for ses = 1:nofsessions;
                if TRL == 7; % Minor adjusment for unequal sessionlengths... 
                   sessions = [sessions, floor(sessionlength * (combinations(comb,ses) -1))+ 1 : floor(sessionlength * combinations(comb,ses))];
                else
                    sessions = [sessions, (combinations(comb,ses)-1)*sessionlength + 1 : combinations(comb,ses)*sessionlength];

                end
            end
            
            RXX = squeeze(mean(covar.Rxx_all(S,sessions,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,sessions,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
            
            others = 1:size(covar.Rxy_at_all,2);
            others(sessions) = [];
            
            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\TRL 7, 8')
            save(['Decoder_Subj' num2str(S) '_Sessions' mat2str(combinations(comb,:)) '_TRL' num2str(TRL) '.mat'], 'avg_dec', 'others')
            
        end
    end
    
end
end

%% TEST THE DECODERS % Same triallength for making the decoder as for evaluating

clear all 

% Adjustable parameters
%numofsessions = 4; % Number of sessions you want the decoder to be built on. 
nofsubj = 4;

% Testing the decoders

for nofsessions = 1:4
    combinations = nchoosek(1:8, nofsessions);
    nofcomb = size(combinations, 1); % Total number of combinations
for TRL = 1:6
    cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate 500')
    load(['TRL' num2str(TRL) '_ISR8_NSR3.mat'])
    noftrials = size(covar.Rxy_at_all,2);
    CA{TRL} = zeros(1,nofcomb,noftrials); CUA{TRL} = zeros(1,nofcomb,noftrials);
    

    for S = 1:nofsubj
        for dec = 1: nofcomb
            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\1, 2, 3, 4')
            load(['Decoder_Subj' num2str(S) '_Sessions' mat2str(combinations(dec,:)) '_TRL' num2str(TRL) '.mat'])
            
            nofothers = length(others);
            
            for HH = 1:nofothers
                recenv = EE{S}(:,:,others(HH))*avg_dec;
                
                CA{TRL}(S, dec, others(HH)) = corr2(recenv,Attended{S}(others(HH),:)');
                CUA{TRL}(S, dec, others(HH)) = corr2(recenv,Unattended{S}(others(HH),:)');
            end
            
            Accuracy(TRL, S, dec) =  (sum( (CA{TRL}(S,dec,:)) - (CUA{TRL}(S,dec,:)) > 0)/(nofothers))*100
        end
    end
end

cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\Results 1, 2, 3, 4')
save(['Results_' num2str(nofsessions)], 'Accuracy', 'CA', 'CUA')
end

%% TEST THE DECODERS % Different triallength for making the decoder as for evaluating
% This script will test all decoders on the triallengths of 5, 10, 15 and
% 20s.

clear all 

% Adjustable parameters

nofsubj = 4;

% Testing the decoders
for nofsessions = 1:4
    combinations = nchoosek(1:8, nofsessions);
    nofcomb = size(combinations, 1); % Total number of combinations    
    
    CA = {}; CUA = {};
    for TRLses = [7, 1, 8, 2] % Komt overeen met 5, 10, 15, 20.
        cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate 500')
        load(['TRL' num2str(TRLses) '_ISR8_NSR3.mat'])
    
        noftrials = size(covar.Rxy_at_all,2);
        sessionlength = size(covar.Rxy_at_all,2)/8
        
        CA{TRLses} = zeros(1,nofcomb,noftrials,8); CUA{TRLses} = zeros(1,nofcomb,noftrials,8);

         for TRLdec = 1:8
            
            for S = 1:nofsubj
                
                for dec = 1: nofcomb
                    cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\1, 2, 3, 4')
                    load(['Decoder_Subj' num2str(S) '_Sessions' mat2str(combinations(dec,:)) '_TRL' num2str(TRLdec) '.mat'])
                    others = [];
                    sessions = [];
                    for ses = 1:nofsessions;
                        if TRLses == 7; % Minor adjusment for unequal sessionlengths... 
                            sessions = [sessions, floor(sessionlength * (combinations(dec,ses) -1))+ 1 : floor(sessionlength * combinations(dec,ses))];
                        else
                            sessions = [sessions, (combinations(dec,ses)-1)*sessionlength + 1 : combinations(dec,ses)*sessionlength];

                        end
                    end
                    
                    others = 1:size(covar.Rxy_at_all,2);
                    others(sessions) = [];
                    
                    nofothers = length(others);
                    
                    for HH = 1:nofothers
                        recenv = EE{S}(:,:,others(HH))*avg_dec;
                
                        CA{TRLses}(S, dec, others(HH), TRLdec) = corr2(recenv,Attended{S}(others(HH),:)');
                        CUA{TRLses}(S, dec, others(HH), TRLdec) = corr2(recenv,Unattended{S}(others(HH),:)');
                    end
                    
                    Accuracy(TRLses, S, dec, TRLdec) =  (sum( (CA{TRLses}(S,dec,:, TRLdec)) - (CUA{TRLses}(S,dec,:, TRLdec)) > 0)/(nofothers))*100
                end
            end
        end
    end
    cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\Decoders\different trial lengths for evaluation and decoders\better')
    save(['Results_' num2str(nofsessions)], 'Accuracy', 'CA', 'CUA')
end
