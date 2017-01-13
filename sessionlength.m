% Calculate session lengths

for TRL = 1:8
    sessionlength = [];
    for S = 1:4
        ses = 0;
        for LR = 1:2
            if LR == 1
                left = true;
            else
                left = false;
            end
            for session = 1:4
                ses = ses +1;
                if left
                   EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\S' num2str(S) '_Pilot_' num2str(session) '_cheap_L.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                else
                   EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\S' num2str(S) '_Pilot_' num2str(session) '_cheap_R.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                end
                [pks,locs] = findpeaks(double(EEG.data(28,:))); % Locate Triggers
                if session == 4 % Auto add trigger after 6min30s because incongruent length of audio stories so missing trigger.
                    locs(2) = locs(1) + ((395)*EEGsampleRate);
                else
                    locs(2) = locs(1) + ((size(data,1)/fs)*EEGsampleRate);
                end
                EEGaad = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
                for ch = 1:24 % EEG
                    [EEGaad(ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, EEGaad(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                    New(ch,:) = resample(double(EEGaad(ch,:)),NSR,EEGsampleRate);
                end
                EEGaad = New;
                EEGsampleRate = NSR;
                % CREATE TRIALS
                window = triallength*NSR; % seconds * Sampling-rate
                [trs] = mod(size(EEGaad,2),window);
                trn = (size(EEGaad,2)-trs)/window; % number of trials possible.
                sessionlength(S, ses) = trn; 
            end
           
        end
    end
         cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\intermediate\Intermediate 500\sessionlengths')
         save(['Sessionlength_TRL' num2str(TRL)], 'sessionlength')

end