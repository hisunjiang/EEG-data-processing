% ----------------------------------
% Data analysis for the benchmark dataset
% Time domain: Bandpass filt -> Trial average -> subject average
% Frequency domain: Bandpass filt -> Trial average -> subject average -> Data tuncation 
%                                   -> FFT
% Sun, Qiang 10/03/2020
% hisunjiang@outlook.com
% Chongqing University EE BCI lab
% ----------------------------------
clear all
close all
clc

Fs = 250;
Wp = [7/Fs*2, 70/Fs*2];
Ws = [4/Fs*2, 75/Fs*2];
[N, Wn]=cheb1ord(Wp, Ws, 0.5, 40);
[b, a] = cheby1(N, 0.5, Wn);
[h, w] = freqz(b,a);
f=(0:length(w)-1)*Fs/2/length(w);
figure
plot(f,abs(h)), grid;
xlabel('Frequency (Hz)')
ylabel('Magnitude')

for Nsubject=1:35
    string=['../data/S', num2str(Nsubject), '.mat'];
    filename1 = string;
    filename2 = '../data/Freq_Phase.mat';
    load (filename1)
    load (filename2)
    % Oz channel, stimulus is 15Hz
    % (1) Trial average
    data = mean(squeeze(data(62,:,8,:)), 2);

    N=2^nextpow2(length(data));  
    T=1/Fs;
    f1=(0:(N-1))/(N*T);  
    % (2) Bandpass filt
    data1 = filtfilt(b,a,data);
    % (3) fft
    amplitude = abs(fft(data,N));
    amplitude = amplitude*2/N;
    amplitude(1) = amplitude(1)/2;
    data2 = amplitude;
        
    AllSubjectDataFilted(:, Nsubject) = data1;
    AllSubjectDataRaw(:, Nsubject) =data2;
end
% (4) Subject average

%% Time domain Analysis (1)+(2)+(4) or (1)+(4)+(2)
AverDataFilted = mean(AllSubjectDataFilted, 2);
figure
subplot(311)
plot(AverDataFilted,'linewidth', 1), xlabel('Time(s)'), ylabel('Amplitude(¦ÌV)')
set(gca, 'Xtick', [0; 400; 1000; 1500])
set(gca,'XTickLabel',{'0','2','4','6'});
line([125 125], [-3 3], 'linewidth', 1, 'color', 'r','linestyle','--')
line([1375 1375], [-3 3], 'linewidth', 1, 'color', 'r','linestyle','--')

%% Frequency domain Analysis (1)+(3)+(4), CAUTION!!! (1)+(4)+(3) is wrong
AverDataRaw = mean(AllSubjectDataRaw, 2);
subplot(312)
plot(f1, AverDataRaw,'linewidth', 1), xlim([5 80]), ylim([0 3])
xlabel('Frequency(Hz)'), ylabel('Amplitude(¦ÌV)')

%% SNR. 1Hz around central frequency, namely, 2048/250=8 points
[~, index] = min(abs(f1-70));
amplitude = AverDataRaw;
SNR=zeros(1, index);
for Nf = 1:index
    if Nf==1
        coef = amplitude(Nf)/mean(amplitude(Nf:Nf+8));
    elseif Nf>1 && Nf<=8
        coef = amplitude(Nf)/mean(amplitude(1:Nf+8));
    elseif Nf>8 && Nf<index-8
        coef = amplitude(Nf)/mean(amplitude(Nf-8:Nf+8));
    elseif Nf>=index-8 && Nf < index
        coef = amplitude(Nf)/mean(amplitude(Nf-8:index));
    else
         coef = amplitude(Nf)/mean(amplitude(Nf-8:Nf)); 
    end
    SNR(Nf) = 20*log10(coef);
end
subplot(313)
plot(f1(1:index), SNR,'linewidth', 1), xlim([5 70])
xlabel('Frequency(Hz)'), ylabel('SNR(dB)')

%% Time & Frequency Analysis
figure
spectrogram(AverDataFilted,hamming(128),127,128,250, 'yaxis'), ylim([0 80]);

