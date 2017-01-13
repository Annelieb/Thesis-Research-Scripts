% {Display boxplots for two different groups, 
for S = 1:12
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr2\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr2.mat'], 'session_Accuracy')
    SA_mimick(S,:) = [squeeze(session_Accuracy(1,:,1)), squeeze(session_Accuracy(1,:,2))];
    load(['C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\onr1\Accuracies_S ' num2str(S) '_NSR20_HP8_LP1_onr1.mat'], 'session_Accuracy')
    SA_onr1(S,:) = [squeeze(session_Accuracy(1,:,1)), squeeze(session_Accuracy(1,:,2))];  
end

load('C:\Users\Annelies\Documents\Studie\Thesis\Final experiments\Analysis\Online experimentes accuracies put together.mat')

f=figure;   
% Boxplot Online 
position_Online = 1:1:12;  
% Define position for 12 Month_O boxplots  
box_Online = boxplot(Session_Accuracy','colors','k','positions',position_Online,'width',0.18); 
set(gca,'XTickLabel',{' '})  % Erase xlabels   
hold on  % Keep the Month_O boxplots on figure overlap the Month_S boxplots   
% Boxplot for the simulated temperature from January to December 
position_Mimick = 1.3:1:12.3;  % Define position for 12 Month_S boxplots  
box_Mimick = boxplot(SA_mimick','colors','b','positions',position_Mimick,'width',0.18);   

position_Onr1 = 1.6:1:12.6;  % Define position for 12 Month_S boxplots  
box_Mimick = boxplot(SA_onr1','colors','r','positions',position_Onr1,'width',0.18);   

hold off   % Insert texts and labels 
 ylabel('Decoding Accuracy [%]') 
 xlabel('Subject')
 
 
% text('Position',[1.1,0],'String','Subject 1') 
% text('Position',[2.1,0],'String','Subject 2') 
% text('Position',[3.1,0],'String','Subject 3') 
% text('Position',[4.1,0],'String','Subject 4') 
% text('Position',[5.1,0],'String','Subject 5') 
% text('Position',[6.1,0],'String','Subject 6') 
% text('Position',[7.1,0],'String','Subject 7') 
% text('Position',[8.1,0],'String','Subject 8') 
% text('Position',[9.1,0],'String','Subject 9') 
% text('Position',[10.1,0],'String','Subject 10') 
% text('Position',[11.1,0],'String','Subject 11') 
% text('Position',[12.1,0],'String','Subject 12') 
% set(gca,'XTickLabel',{''});   % To hide outliers 
% out_O = box_O(end,~isnan(box_O(end,:)));  
% delete(out_O)  
% out_S = box_S(end,~isnan(box_S(end,:)));  
% delete(out_S)
% legend('Online', 'Mimick', 'Offline')