function [num] = Rbp(L, H, order, fs, Signal)
% Bandpass Butter Filter Rob Zink 2015
Signal = double(Signal);

[B,A] = butter(order,H/(fs/2),'low');
Signal = filtfilt(B,A,Signal);

[B,A] = butter(order,L/(fs/2),'high');
Signal = filtfilt(B,A,Signal);

num = Signal;
end