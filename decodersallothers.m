%% Decoder Pilot study
for TRL = 1:6
    cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\Investigate LP and HP values')
    load(['LP1_HP8_TRL' num2str(TRL) '_ISR8_NSR3.mat'], 'covar')
    Rxx = []; Rxy = [];
   
    for S = [1, 3, 4]
        Rxx = cat(1, Rxx, squeeze(covar.Rxx_all(S,:,:,:)));
        Rxy = cat(1, Rxy, squeeze(covar.Rxy_at_all(S,:,:)));
    end
    
    RXX = squeeze(mean(Rxx,1));
    RXY = squeeze(mean(Rxy,1));
    
    avg_dec{TRL} = RXX \ RXY'; % LS solution.
end
%%
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
save(['Pilot_decoder_134'], 'avg_dec')

%% Decoder of all other
cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
subj = 1:12
for S = 1:12
    Rxx = cell(1,6); Rxy = cell(11,6);
    others = subj(2:end)
    for o = others
        cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1')
        load(['Accuracies_S ' num2str(o) '_NSR20_HP8_LP1_onr1.mat'], 'Rxx_all', 'Rxy_at_all')
        for TRL = 1:6
            Rxx{TRL} = cat(1, Rxx{TRL}, Rxx_all{TRL});
            Rxy{TRL} = cat(1, Rxy{TRL}, Rxy_at_all{TRL});
        end
    end
    subj = [subj(2:end), subj(1)]
    
    for TRL = 1:6
        RXX = squeeze(mean(Rxx{TRL}));
        RXY = squeeze(mean(Rxy{TRL}));
        avg_dec{TRL} = RXX\RXY';
    end
    cd('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis')
    save(['Decoder all other S' num2str(S)])
end
    