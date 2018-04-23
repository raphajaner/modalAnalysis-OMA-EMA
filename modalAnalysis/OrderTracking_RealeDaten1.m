%% 1. Clean-Up
clc, clearvars, clf

%% 2. Initialization

% Pfad hinzufügen
addpath('../../Daten_Ordertracking');
addpath('./functions');

% Sensordaten einladen
load Beschleunigungssensoren_2Stueck.mat

% Drehzahl Trommel (?) einladen
% load Trommeldrehzahl_Schnell_Geglaettet.mat
% load Trommeldrehzahl_Schnell
load Trommeldrehzahl_Langsam.mat

% Sampling Frequenz
fs = 2048; 

% Problem: Zeitpunkte der Trommeldrehzahlen sind doppelt vorhanden, damit
% sind unterschiedliche Werte den Zeitpunkten zugeordnet. Die Funktion 
% griddedInterpolant benötigt allerdings unique Datenpunkte. Daher wird nur
% jeder dritte Punkte ausgewählt und die Zwischenpunkte interpoliert.

TrommeldrehzahlS = TrommeldrehzahlS(1:size(TrommeldrehzahlS)-1);
TrommeldrehzahlS_Zeit = TrommeldrehzahlS_Zeit(1:size(TrommeldrehzahlS_Zeit)-1);

% Interpoliertes Element
TrommeldrehzahlS_interpol = griddedInterpolant(TrommeldrehzahlS_Zeit,TrommeldrehzahlS);

% Auswertung des interpolierten Elements an den Zeitpunkten, an denen die
% Sensoren ausgewertet wurde. Umrechnung von [1/s] in [rpm]
Trommeldrehzahl_RMP = TrommeldrehzahlS_interpol(Y_Zeit)*60;
Trommeldrehzahl_RMPSmooth = smoothdata(Trommeldrehzahl_RMP,1);

% Normalerweise hier Aufruf der implementierten Funktion orderTracking.
% Hierbei jedoch das Problem, dass bei pcg der iterative Funktionsaufruf
% nicht konvergiert.
% [mag, freqLine, orderlist] = orderTracking(Y_6807_X,fs,Trommeldrehzahl_RMP,2);

% Probeweise neu aufgebaute Prozedure als Test. Hierbei wird der von Scot
% McNeill bereitgestellte Algorithmus verwendet. Dieser bildet Inverse,
% anstatt pcg zu verwendeten.

% Visualisierung FreqMap
rpmfreqmap(Y_6807_X,fs,Trommeldrehzahl_RMP)
% Bestimmung der gesuchten Ordnungen
[spec,order] = orderspectrum(Y_6807_X,fs,Trommeldrehzahl_RMP);
[~, orderlist] = findpeaks(spec, order,'SortStr','descend','Npeaks', 10);
orderlist = round(orderlist,1);
display(orderlist);

% Vold-Kalman-Filter 
[x_vk2,bw,T,xr_vk2] = vk2(Y_6807_X,Trommeldrehzahl_RMP/60*orderlist(3),fs,1.5*32e6,1);

% Eigenfrequenzen identifizieren
[~, omega] = findpeaks(abs(x_vk2), Trommeldrehzahl_RMP(end:-1:1)/60*orderlist(3),'SortStr','descend','Npeaks', 5);
findpeaks(abs(x_vk2), Trommeldrehzahl_RMP(end:-1:1)/60*orderlist(3),'SortStr','descend','Npeaks', 5);
display(omega);