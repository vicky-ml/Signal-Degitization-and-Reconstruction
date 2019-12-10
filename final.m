%Author : Vignesh Waren Sunder
%Signal Digitisation and Reconstruction
%Coded using MATLAB R2019a- academic use
%GitHub: vicky-ML


close all
clc

%Generating Signal
t = 0:0.001:2; % set the time domain range (1 to 2 with 1000 steps)
f1 = 6; 
f2 = 9;
x = sin(2*pi*f1*t)+sin(2*pi*f2*t); %Sanalog
figure
subplot(211)
plot(t,x);
title("Sanalog")
xlabel("Time(s)");
Fs = 100; %Our Sampling Frequency
Ts = 1/Fs;
n = 0:Ts:2; %set the time domain range
x_sampled= sin(2*pi*f1*n)+sin(2*pi*f2*n); %Ssam
subplot(212)
stem(n,x_sampled,'fill');
title("Ssam");
xlabel("Time(s)  Fs=100Hz");





%Quantisation
maxx= max(x_sampled); %max value of Ssam
minn = min(x_sampled); % min value of Ssam
range = (maxx-minn)/16; % getting range with appropriate quatisation level

%we perform "for" loop to get the y_axis value for each step in x_axis and
%then round it to the nearest integer
for n1 = 1:length(x_sampled)
    
    x_quant=round((x_sampled-minn)/range)*range+minn;
    
end
%Signal to quantisation noise ratio
powerinsig = x_sampled.^2; %Power in Signal
powerinnoise = (x_sampled-x_quant).^2; %Power in Noise
sqnr = 10*log10 (powerinsig/powerinnoise); %SQNR

figure
plot(t,x);
hold
stairs(n,x_quant)
title('Quantised')
xlabel("Time(s)");

%fft
%The value of Fs changes. Fs=100 for FFT of Ssamp | Fs=1000 for FFT of
%Sanalog. Kindly change the value when running the code.
xft = fft(x_sampled(1:length(x_sampled)-1)); %perform fft with a point delay(-1)
arr=0:Fs/(2*Fs):Fs-Fs/(2*Fs); %To get the range of the frequency spectrum
%(arr is 0 to 200 with 2000 steps) gives for one cycle
%We generated and plotted Frequency spectrum only for One Cycle so that we  
%can easily perform FFTShift and extract the required frequency component
xabs=abs(xft); %gives the magnitude
figure
subplot(211)
plot(arr,xabs/Fs)
title('FFT of Sanalog')
xlabel("Frequency(Hz)");
ylabel("Magnitude");
xxftshift = fftshift(xft); %perform FFTShift
subplot(212)
plot(arr-0.5*Fs, abs(xxftshift)/Fs)%plot fftshift with desired range
title("FFTShift of Sanalog");
xlabel("Frequency(Hz)");
ylabel("Magnitude");

%Low Pass Filter
lpfilter = zeros (1,length(xxftshift)); %zero magnitude for this range

for number1 = 1:length(xxftshift) %magnitude is 1 for these range of points
    
    if number1< 0.5*length(xxftshift)
        lpfilter(number1)=1;
    end
end
figure
plot(arr,xabs/Fs.*lpfilter) %perform cyclic convolution by using point wise
%multiplicaiton in Frequency domain. The required frequency signal is 
%extracted
title("filtered spectrum Sanalog");
xlabel("Frequency(Hz)");
ylabel("Magnitude");

figure
iff = ifft(xft/Fs.*lpfilter*2); %perform ifft with extracted signal from
%the lpfilter
%iff = ifft(xft/Fs);
plot(n(1:length(n)-1),iff*Fs, 'r') %Plot the reconstructed Signal
hold
%stem(n,x_sampled,'fill'); %to get a stemed reconstructed signal output
%hold 
plot(t,x, 'k')
title("Reconstructed Signal and Sanalog")
xlabel("Time (s)");
ylabel("Amplitude");
