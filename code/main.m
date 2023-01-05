clc
close all
clear all

%% Loading ECG and pulse oximetry (PO) signals

signals = load('samples/ECGNormal+PO.txt');
ecg = signals(:,1); % ECG sample
po1 = signals(:,2); % PO sample for wavelength 600..700 nm
po2 = signals(:,3); % PO sample for wavelength 810..960 nm
n = length(ecg);
fs = 250;  % Sampling frequency equals 250 Hz
T = 1/fs;
tMax = (n-1)/fs;
t = 0:T:tMax; % Array of timeline

clear signals

%% Filter ECG

fNet = 60; % Network tip-off rate = 60 Hz
ecgToFilter = ecg + 200*sin(2*pi*fNet*t)'; % Adding network noise

load('filters/pitch60.mat'); % Loading pitch filter coeffs (60 Hz)
ecgFiltered = filter(Num,Den,ecgToFilter); % Filter ecg

figure(1) % Building plot

subplot(2,1,1)
plot(t,ecgToFilter)
title("Before filtering")
xlabel("Time, s")
ylabel("Magnitude, mV")
grid on
ylim([-300 600])

subplot(2,1,2)
plot(t,ecgFiltered)
title("After filtering")
xlabel("Time, s")
ylabel("Magnitude, mV")
grid on
ylim([-300 600])

clear Num Den fNet ecgFiltered ecgToFilter;

%% Spectrum analysis

df = fs/n;
f = 0:df:fs-df; % Array of spectrum frequencies
ecgSpectrumPower = abs(fft(ecg)); % Spectral power density
ecgSpectrumPower = ecgSpectrumPower(1:n/2+1);
ecgSpectrum = ecgSpectrumPower / n; % Amplitude spectrum
ecgSpectrum(2:end) = 2*ecgSpectrum(2:end);

figure(2) % Building spectrum plot

subplot(2,1,1)
plot(f(1:n/2+1),ecgSpectrumPower)
title("Spectral power density")
xlabel("Frequency, Hz")
ylabel("Power")
grid on
xlim([0 fs/2])

subplot(2,1,2)
plot(f(1:n/2+1),ecgSpectrum)
title("Amplitude spectrum")
xlabel("Frequency, Hz")
ylabel("Amplitude, mV")
grid on
xlim([0 fs/2])

clear ecgSpectrum ecgSpectrumPower f df;

%% Calculation of arterial blood oxygen saturation

ecgDiff = diff(ecg); % Derivative of the ECG signal
absEcgDiff = abs(ecgDiff); % Absolute derivative of the ECG signal

threshold = max(absEcgDiff) / 4; % Threshold for QRS detect algorithm
refT = round(0.25 / T); % Refractory period for QRS detect algorithm

% QRS detect algorithm
j = 0;
k = 0;
for i = 1:n-1
    j = j+1;
    if (absEcgDiff(i) > threshold)&&(j > refT)
            k = k+1;
            QRS(k) = i;
            j = 0;
    end
end

% Calculating average values of the range of oscillation of PO signals
d1 = 0;
d2 = 0;
for i = 1:(k-1)
    m1 = po1(QRS(i):QRS(i+1));
    d1 = d1 + (max(m1)-min(m1));
    m2 = po2(QRS(i):QRS(i+1));
    d2 = d2 + (max(m2)-min(m2));
end
d1 = d1 / (k - 1);
d2 = d2 / (k - 1);

alpha = d2 / d1;
satO2 = (0.872 - 0.16 * alpha) * 100 / (0.14 * alpha + 0.754); % Saturation

figure(3) % Building plot

subplot(4,1,1)
plot(t, ecg)
title("ECG")
xlabel("time, s")
ylabel("Amplitude, mV")
grid on

subplot(4,1,2)
plot(t(2:end), ecgDiff)
title("Derivative ECG")
xlabel("time, s")
grid on

subplot(4,1,3)
plot(t(2:end), absEcgDiff)
hold on
plot([0 tMax], [threshold threshold])
title("Absolute derivative ECG")
xlabel("time, s")
grid on
subplot(4,1,4)
plot(t, po1)
hold on
plot(t, po2)
for i=1:k
    line([t(QRS(i)) t(QRS(i))], [min([po1; po2]) max([po1; po2])],'LineStyle','--')
    scatter(t(QRS(i)),po1(QRS(i)),'blue')
    scatter(t(QRS(i)),po2(QRS(i)),'red')
end
title("Pulse oximetry")
xlabel("time, s")
ylabel("Amplitude, mV")
grid on

disp(['Saturation: ', num2str(satO2)]); % Display saturation value in CW

clear absEcgDiff alpha d1 d2 ecgDiff i j k m1 m2 QRS refT satO2 threshold;

%% Сorrelation analysis

ecgCorr = xcorr(ecg, 'coeff', n/2); % EСG autocorrelation sequence 
poCorr = xcorr(po1, po2, 'coeff', n/2); % PO cross-correlation
tCorr =-tMax/2:T:tMax/2+T; 

figure(4) % Building plot

subplot(2,1,1)
plot(tCorr, ecgCorr)
title("EСG autocorrelation sequence")
xlabel("time, s")
xlim([-tMax/2 tMax/2])
grid on

subplot(2,1,2)
plot(tCorr, poCorr)
title("PO cross-correlation")
xlabel("time, s")
xlim([-tMax/2 tMax/2])
grid on

clear ecgCorr poCorr tCorr;
clear ecg n po1 po2 t tMax;

%% Classification of ECG with сorrelation analysis

ecg = load('samples/ECGPathology.txt');
markers = load('samples/ECGPathologyMarking.txt'); % Index and marker
n = length(ecg);
tMax = (n-1)/fs;
t = 0:T:tMax; % Array of timeline
nQRS = length(markers);

% Loading QRS segments for 2 norms and pathologies for visualisation
n0=1;
p0=2;
n1=3;
p1=8;
l=51;
SMax=(l-1)/2;
tc=-T*SMax:T:T*SMax;

ind0=markers(n0,1);
QRS_n0=detrend(ecg(ind0-SMax:ind0+SMax));

ind1=markers(n1,1);
QRS_n1=detrend(ecg(ind1-SMax:ind1+SMax));

ipd0=markers(p0,1);
QRS_p0=detrend(ecg(ipd0-SMax:ipd0+SMax));

ipd1=markers(p1,1);
QRS_p1=detrend(ecg(ipd1-SMax:ipd1+SMax));

% Calculating correlation values
for i=1:nQRS
    QRS=ecg(markers(i,1)-SMax:markers(i,1)+SMax);
    QRS=detrend(QRS);
    cn=xcorr(QRS_n0,QRS,5,'coeff');
    cp=xcorr(QRS_p0,QRS,5,'coeff');
    cnMax(i)=max(cn);
    cpMax(i)=max(cp);
end

% QRS classification using correlation values 
corrLimit = 0.95; % Correlation limit for algorithm accuracy
for i=1:nQRS
    if (cnMax(i)>=corrLimit)&&(cpMax(i)<corrLimit)
        classes(i)=1; % Classes created by algorithm
    elseif (cpMax(i)>=corrLimit)&&(cnMax(i)<corrLimit)
        classes(i)=2;
    else
        classes(i)=0;
    end
end

figure(5) % Building plot

subplot(4,3,1:3)
plot(t, ecg)
title("ECG")
xlabel("Time, s")
ylabel("Magnitude, V")
grid on

subplot(4,3,4:6)
for i=1:nQRS
    x=markers(i,1)*T;
    if (markers(i,2) == 1) 
        marker = 'N';
    else
        marker = 'P';
    end
    if (classes(i) == 1) 
        class = 'N';
    else
        class = 'P';
    end
    text(x,0.7,marker,'Color','b');
    text(x,0.3,class,'Color','r');
end
title("QRS markers")
xlabel("Time, s")
xlim([0 tMax])
ylim([0 1])
grid on

subplot(4,3,7)
plot(tc,QRS_n0);
title("Norm")
xlabel("Time, s")
ylabel("Magnitude, V")
ylim([-0.5 1])
grid on

subplot(4,3,8)
plot(tc,QRS_n0,'blue')
hold on
plot(tc,QRS_n1,'red')
plot(tc,xcorr(QRS_n0,QRS_n1,SMax,'coeff'),':black', 'LineWidth', 2)
legend('Norm', 'Norm', 'Correlation')
ylim([-1 1])
grid on
title('N-N')

subplot(4,3,9)
plot(tc,QRS_n0,'blue')
hold on
plot(tc,QRS_p1,'red')
plot(tc,xcorr(QRS_n0,QRS_p1,SMax,'coeff'),':black', 'LineWidth', 2)
legend('Norm', 'Pathology', 'Correlation')
ylim([-1 1])
grid on
title('N-P')

subplot(4,3,10)
plot(tc,QRS_p0);
title("Pathology")
xlabel("Time, s")
ylabel("Magnitude, V")
ylim([-0.5 1])
grid on

subplot(4,3,11)
plot(tc,QRS_p0,'blue')
hold on
plot(tc,QRS_n1,'red')
plot(tc,xcorr(QRS_p0,QRS_n1,SMax,'coeff'),':black', 'LineWidth', 2)
legend('Pathology', 'Norm', 'Correlation')
ylim([-1 1])
grid on
title('P-N')

subplot(4,3,12)
plot(tc,QRS_p0,'blue')
hold on
plot(tc,QRS_p1,'red')
plot(tc,xcorr(QRS_p0,QRS_p1,SMax,'coeff'),':black', 'LineWidth', 2)
legend('Pathology', 'Pathology', 'Correlation')
ylim([-1 1])
grid on
title('P-P')

clear class classes cn cnMax corrLimit cp cpMax i ind0 ind1 ipd0 ipd0 ipd1
clear l marker markers n0 n1 nQRS p0 p1 QRS QRS_n0 QRS_n1 QRS_p0 QRS_p1
clear SMax tc x