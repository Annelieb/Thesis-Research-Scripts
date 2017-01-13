% Comparing the bp filters for S4 TRL 6
LP = 2; HP = 9; intermediatSamplingRate = 500; triallength = 60; 
S = 4; ses = 0; power = 0.6;

for LR = 1:2 % Left - Right loop
        if LR == 1
            left = true;
        else
            left = false;
        end

        for session = 1:4 % 4 right ear session - Train on each one separate
            ses = ses + 1;
            data = [];
            [data2,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track1_dry.wav']);
            [data3,fs] = audioread(['C:\Users\Annelies\Documents\Studie\Thesis\stimuli\part' num2str(session) '_track2_dry.wav']);
            data = cat(2,data2(:,1),data3(:,1));
            
                    if left
                        EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\off line pilot\S' num2str(S) '_Pilot_' num2str(session) '_cheap_L.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                    else
                        EEG = pop_loadxdf(['C:\Users\Annelies\Documents\Studie\Thesis\script\raw_data\off line pilot\S' num2str(S) '_Pilot_' num2str(session) '_cheap_R.xdf'], 'streamname', 'EEG', 'streamtype', 'EEG', 'exclude_markerstreams', {});
                    end
            
            EEGsampleRate = 500;
            
            % DATA TRUNCATION & DOWNSAMPLING
            [pks,locs] = findpeaks(double(EEG.data(28,:))); % Locate Triggers
 
            if session == 4 % Auto add trigger after 6min30s because incongruent length of audio stories so missing trigger.
                locs(2) = locs(1) + ((395)*EEGsampleRate);
            else
                locs(2) = locs(1) + ((size(data,1)/fs)*EEGsampleRate);
            end
            EEGaad{ses} = EEG.data(1:24,locs(1):locs(2)); % Truncate EEG
            
            AudioL = data(:,1); AudioR = data(:,2);
            data = [];
            
            % Rob resamples before making the envelope
                AudioL = resample(AudioL,EEGsampleRate,fs); % Downsampled audio per ear
                AudioR = resample(AudioR,EEGsampleRate,fs);
                fs = EEGsampleRate;  
            
            % ENVELOPES
            % Powerlaw
                envelopeL{ses} = abs(AudioL).^power;
                envelopeR{ses} = abs(AudioR).^power;
                        
            % FILTERING
            for ch = 1:24 % EEG
                [EEGaad_old{ses}(ch,:)] = Rbp_old(LP, HP, 6, EEGsampleRate, EEGaad{ses}(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
                [EEGaad_new{ses}(ch,:)] = Rbp(LP, HP, 6, EEGsampleRate, EEGaad{ses}(ch,:)); % Band pass % OPM: In de nieuwe bp filter is order niet meer van belang! 
            end

            [envelopeL_old{ses}] = Rbp_old(LP, HP, 6, fs, envelopeL{ses}(:,1)); % Band pass
            [envelopeR_old{ses}] = Rbp_old(LP, HP, 6, fs, envelopeR{ses}(:,1)); % Band pass
            
            [envelopeL_new{ses}] = Rbp(LP, HP, 6, fs, envelopeL{ses}(:,1)); % Band pass
            [envelopeR_new{ses}] = Rbp(LP, HP, 6, fs, envelopeR{ses}(:,1)); % Band pass
        end
end

%%
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter')
save('CompareFilters', 'envelopeL', 'envelopeR', 'envelopeL_old', 'envelopeR_old', 'envelopeL_new', 'envelopeR_new', 'EEGaad', 'EEGaad_old', 'EEGaad_new')
%% Make plots
for ses = 1
    % audio
    figure
    plot(envelopeL{ses}(1:3000))
    hold on
    plot(envelopeL_old{ses}(1:3000))
    plot(envelopeL_new{ses}(1:3000))
    hold off
    legend('original', 'old', 'new')
end