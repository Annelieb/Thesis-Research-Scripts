%% Make decoders
% Initial parameters

close all
clear all

nofsubj = 3;
lag = [100 400];
LP = 2;
HP = 9;
NSR = 20;
triallength = 10

envelopemethod = 'powerlaw'; %or hilbert
power = 0.6; % Why? Ask Neetha... From literature?

%% Subject accuracy
cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []')
load(['TRL1_ISR8_NSR3.mat'])
Subject_Accuracy = [];
    
    for S = 1:4
        
            RXX = squeeze(mean(covar.Rxx_all(S,:,:,:),2)); % plain averaging of cov matrices
            RXY_ATT = squeeze(mean(covar.Rxy_at_all(S,:,:),2));
        
            avg_dec = RXX \ RXY_ATT; % LS solution.
     
            cd('C:\Users\Annelies\Documents\Studie\Thesis\script\All_data\New bp filter\audio HL = []\General decoders')

     save(['TRL1_Decoder_S' num2str(S) '.mat'], 'avg_dec')
    end

       
