function [mag, freqLine, orderlist] = orderTracking(vibData, fs, rpmVec, numOrders)
% Bestimmung der Ordnungen und deren rms-Wert
[spec,order] = orderspectrum(vibData,fs,rpmVec);

% 1. Subplot. Die betragsmäßig größten Ordnungen werden identifiziert und
% visualisiert (numOrders bestimmt Anzahl).
subplot(1,2,1)
findpeaks(spec, order,'SortStr','descend','Npeaks', numOrders);

% Bestimmung der gesuchten Ordnungen
[~, orderlist] = findpeaks(spec, order,'SortStr','descend','Npeaks', numOrders);

% Spaltenindizees für ordertrack
rpmrefidx = linspace(1, numOrders, numOrders);

% Die am meisten zum System beitragenden Ordnungen werden getrackt und
% dazugehörend magnitude, rpm und time bestimmt.
[mag, rpm, ~] = ordertrack(vibData, fs, [rpmVec rpmVec], orderlist, rpmrefidx, 'Decouple', true);

% Transponieren der Matrizen/Vektoren
mag = mag';
orderlist = orderlist';
rpm = rpm';

% Frequenzband der jeweiligen Ordnungen
freqLine = rpm/60.*orderlist;

% 2. Subplot. Visualisiert das Frequenzspektrum
subplot(1,2,2)
hold on

% Die drei größten Peaks werden indentifiziert
for numOrders_run = 1 : numOrders
    findpeaks(mag(: ,numOrders_run), freqLine(:, numOrders_run), 'SortStr','descend','Npeaks', 3)
end

end

