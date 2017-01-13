% Overnight 17/07

% Channel reduction
% average over all subjects to get an overal removedchannels
% Get the covariances
subjects = 1:12;
startleft = [1 0 1 0 1 0 0 1 0 1 0 1];

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6;

lag1 = 100;
lag2 = 400;                                                    
lp = 1; % in pilot 2 tot 9 Hz bekeken
hp = 8;
nsr = 20
trll = [10 30 60] % 20 30 40 50 60]; % Eventueel ook voor andere trll? 
intSamp = 500;

l1 = lag1; l2 = lag2; LP = lp; HP = hp; NSR = nsr; intermediateSampleRate = intSamp;

nofreallags = 7; % Depends on NSR and chosen upper and lower limit of lags
for TRL = 1:length(trll)
    
    %Channel reduction
    remainingchannels = 1:24;

    for reduction = 1:24
        
        nofremainingch = length(remainingchannels); 
        
        if isequal(nofremainingch, 24-reduction+1) == 0
            error('Mistake occured')
        end
        
        for S = subjects
            EE = {}; Attended = {}; Unattended = {}; Rxx_all ={}; Rxy_at_all = {};
            Elags = {};
    
            cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
            load(['Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'],'EE', 'Attended', 'Unattended', 'Rxx_all')
            
            for i = 1:size(Rxx_all,2)
                noftrl(i) = size(Rxx_all{i},1);
            end
            clear Rxx_all
            
            window = size(Attended{TRL,1},1)
            
            % Covariance Matrices
            noftrials = noftrl(TRL);
            
            for trial = 1:noftrials
                X = [];
                Elags{TRL,trial} = reshape(EE{TRL,trial},[window,nofreallags,24]); % You can use this for the seperate lags; or for channel reduction!
                X = reshape(Elags{TRL,trial}(:,:,remainingchannels), [window, nofreallags*nofremainingch]);
                
                EEred{TRL,trial, reduction} = X;
                
                nofsamples = size(X,1);
                Rxx_all{TRL,reduction}(trial,:,:) = (X'*X)/nofsamples;

                Rxy_at_all{TRL, reduction}(trial,:) = (X'*Attended{TRL,trial})/nofsamples;
            end

            %Subject Accuracy
            clear CA CUA train
            train = (1:noftrials);
            for HH = 1:noftrials % repeat for all trials - leave one out
                clear RXX RXY_ATT
                trials2 = train(1:end-1);

                RXX = squeeze(mean(Rxx_all{TRL, reduction}(trials2,:,:),1)); % plain averaging of cov matrices
                RXY_ATT = squeeze(mean(Rxy_at_all{TRL, reduction}(trials2,:),1));

                avg_dec{TRL, reduction}(train(end),:) = RXX \ RXY_ATT'; % LS solution.
                recenv{TRL, train(end), reduction} = EEred{TRL, train(end), reduction} * squeeze(avg_dec{TRL, reduction}(train(end),:))';

                CA(train(end)) = corr2(recenv{TRL, train(end), reduction},Attended{TRL,train(end)}); % Attended
                CUA(train(end)) = corr2(recenv{TRL, train(end), reduction},Unattended{TRL, train(end)}); % Unattended
               
                
                % For the channel reduction the decoder needs to be
                % rescaled:
                I = inv(RXX);
                D = diag(I);
                F{TRL, reduction}(S,train(end),:) = (avg_dec{TRL, reduction}(train(end),:).^2)./D'; 
                
                train = [train(end) train(1:end-1)];
            end
            
            subject_CA{S, TRL,reduction} = CA;
            subject_CUA{S, TRL,reduction} = CUA;
            subject_difference{S, TRL, reduction} = CA - CUA;
            subject_Accuracy{TRL}(S, reduction) = (sum(CA-CUA> 0)/noftrials)*100 
            
            clear EE EEred Attended Unattended Rxx_all Rxy_at_all subject_CA subject_CUA subject_difference avg_dec recenv 
        end
        
        FT = [];
        for S = 1:12
            FT = cat(1,FT, squeeze(F{TRL, reduction}(S,:,:)));
        end
        
        % Channel reduction: decide which channel should be eliminated
        % 1) Average the leave one out decoder
        Mdec{TRL, reduction} = squeeze(mean(FT,1));
        % 2) Reshape -> 7x24 matrix
        RMdec{TRL, reduction} = reshape(Mdec{TRL, reduction},[nofreallags, nofremainingch]);
        % 3) Mean over lags
        MRMdec{TRL, reduction} = mean(abs(RMdec{TRL, reduction}), 1);
        % 5) Find channel weight closest to 0 and remove it
        [M, I] = min(MRMdec{TRL, reduction})
        removedchannels{TRL}(reduction) = remainingchannels(I);
        remainingchannels(I) = [];

    end
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis'])
    save(['CR_allS_trl_' num2str(TRL) '.mat'], 'subject_Accuracy', 'F', 'Mdec', 'RMdec', 'MRMdec', 'removedchannels')
    clear subject_Accuracy F Mdec RMdec MRMdec removedchannels
end

% Initial parameters

subjects = 1:12;

% Get the covariances
envelopemethod = 'powerlaw'; %or hilbert
power = 0.6;

lag1 = 100;
lag2 = 400;                                                    
lp = 1; % in pilot 2 tot 9 Hz bekeken
hp = 8;
nsr = 20
trll = [30 60] % 20 30 40 50 60]; % Eventueel ook voor andere trll? 
intSamp = 500;

l1 = lag1; l2 = lag2; LP = lp; HP = hp; NSR = nsr; intermediateSampleRate = intSamp;

nofreallags = 7; % Depends on NSR and chosen upper and lower limit of lags
for S = subjects
    EE = {}; Attended = {}; Unattended = {}; Rxx_all ={}; Rxy_at_all = {};
    Elags = {};
    
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
    load(['Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'],'EE', 'Attended', 'Unattended', 'Rxx_all')
    for i = 1:size(Rxx_all,2)
        noftrl(i) = size(Rxx_all{i},1);
    end
    clear Rxx_all
    for TRL = 1:length(trll)
        window = size(Attended{TRL,1},1)

        %Channel reduction
        remainingchannels = 1:24;
        for reduction = 1:24 % number of removed channels = reduction - 1
            nofremainingch = length(remainingchannels); 
            if isequal(nofremainingch, 24-reduction+1) == 0
                error('Mistake occured')
            end
            % Covariance Matrices
            noftrials = noftrl(TRL);
            for trial = 1:noftrials
                X = [];
                Elags{TRL,trial} = reshape(EE{TRL,trial},[window,nofreallags,24]); % You can use this for the seperate lags; or for channel reduction!
                X = reshape(Elags{TRL,trial}(:,:,remainingchannels), [window, nofreallags*nofremainingch]);
                
                EEred{TRL,trial, reduction} = X;
                
                nofsamples = size(X,1);
                Rxx_all{TRL,reduction}(trial,:,:) = (X'*X)/nofsamples;

                Rxy_at_all{TRL, reduction}(trial,:) = (X'*Attended{TRL,trial})/nofsamples;
            end

            %Subject Accuracy
            clear CA CUA train
            train = (1:noftrials);
            for HH = 1:noftrials % repeat for all trials - leave one out
                clear RXX RXY_ATT
                trials2 = train(1:end-1);

                RXX = squeeze(mean(Rxx_all{TRL, reduction}(trials2,:,:),1)); % plain averaging of cov matrices
                RXY_ATT = squeeze(mean(Rxy_at_all{TRL, reduction}(trials2,:),1));

                avg_dec{TRL, reduction}(train(end),:) = RXX \ RXY_ATT'; % LS solution.
                recenv{TRL, train(end), reduction} = EEred{TRL, train(end), reduction} * squeeze(avg_dec{TRL, reduction}(train(end),:))';

                CA(train(end)) = corr2(recenv{TRL, train(end), reduction},Attended{TRL,train(end)}); % Attended
                CUA(train(end)) = corr2(recenv{TRL, train(end), reduction},Unattended{TRL, train(end)}); % Unattended
               
                
                % For the channel reduction the decoder needs to be
                % rescaled:
                I = inv(RXX);
                D = diag(I);
                F{TRL, reduction}(train(end),:) = (avg_dec{TRL, reduction}(train(end),:).^2)./D'; 
                
                train = [train(end) train(1:end-1)];
            end
            
            subject_CA{TRL,reduction} = CA;
            subject_CUA{TRL,reduction} = CUA;
            subject_difference{TRL, reduction} = CA - CUA;
            subject_Accuracy(TRL, reduction) = (sum(CA-CUA> 0)/noftrials)*100 
            
            % Channel reduction: decide which channel should be eliminated
            % 1) Average the leave one out decoder
            Mdec{TRL, reduction} = mean(F{TRL, reduction},1);
            % 2) Reshape -> 7x24 matrix
            RMdec{TRL, reduction} = reshape(Mdec{TRL, reduction},[nofreallags, nofremainingch]);
            % 3) Mean over lags
            MRMdec{TRL, reduction} = mean(abs(RMdec{TRL, reduction}), 1);
            % 5) Find channel weight closest to 0 and remove it
            [M, I] = min(MRMdec{TRL, reduction})
            removedchannels{TRL}(reduction) = remainingchannels(I);
            remainingchannels(I) = [];
        end
    end
    cd(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis'])
    save(['Channel_reduction2_S ' num2str(S) '_trl3060.mat'], 'EE', 'EEred', 'Rxx_all', 'Rxy_at_all', 'Attended', 'Unattended', 'subject_CA',...
        'subject_CUA', 'subject_difference', 'subject_Accuracy', 'avg_dec', 'recenv', 'Mdec', 'RMdec', 'MRMdec', 'removedchannels')
    clear EE EEred Attended Unattended Rxx_all Rxy_at_all subject_CA subject_CUA subject_difference subject_Accuracy avg_dec recenv Mdec RMdec MRMdec removedchannels
end

