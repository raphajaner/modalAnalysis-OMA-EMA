%% 1. Clean-Up
clc, clearvars, clf

addpath('./functions');

%% 2. Synthetisches Messdaten generieren
% Basisparameter
fs = 300;
nsamp = 3*fs;
time_vec = linspace(0, nsamp, nsamp+1)*1/fs;

% Drehzahlvektor
rpmVec = linspace(10, 100, nsamp)'*60;

% order 1.2 hat 3 Eigenfrequenzen. order 0.8 hat 1 Eigenfrequenz
% (physikalisch ist es nicht möglich, dass orders unterschiedliche
% Eigenfrequenzen haben, zum Test der Funktion jedoch hinreichend.
x_1 = [2 4 3 2]*sqrt(2).*cos(2*pi*cumtrapz([1.2*rpmVec 1.2*rpmVec 1.2*rpmVec 0.8*rpmVec]/60)/fs);

% Teil, der die Eigenfrequenzen darstellt
rs_1 = [1+1./(1+linspace(-10,10,nsamp).^4)'/2 1+1./(1+linspace(-6,14,nsamp).^4)'/2 ...
    1+1./(1+linspace(-2,18,nsamp).^4)'/2 1+1./(1+linspace(-5,15,nsamp).^4)'/2];

% Überlagerung der erzeuften Schwingungen
vibData = sum(rs_1.*x_1,2);

%% 3. Ordertracking Analyse

% Erzeugen der Frequenzmap: 3D-Map mit Frequenz über rpm (entlang
% der order) mit Amplitude (rms)

% rpmfreqmap(vibData,fs,rpmVec,'OverlapPercent',0.1,'Window','hamming')

% Auswertung entlang der größten Ordnung
[mag, freqLine, orderlist] = orderTracking(vibData,fs,rpmVec);

[~, omega] = findpeaks(mag(:,2), freqLine(:,1),'SortStr','descend','Npeaks', 3);

display(omega);







