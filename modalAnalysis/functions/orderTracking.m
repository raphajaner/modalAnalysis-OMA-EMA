function [mag, freqLine, orderlist] = orderTracking(vibData,fs,rpmVec)
% Bestimmung der Ordnungen und deren rms-Wert
[spec,order] = orderspectrum(vibData,fs,rpmVec);

% 1. Subplot. Die vier betragsmäßig größten orders werden identifiziert und
% visualisiert
subplot(1,2,1)
findpeaks(spec, order,'SortStr','descend','Npeaks', 2);

[~, orderlist] = findpeaks(spec, order,'SortStr','descend','Npeaks', 2);

% Die am meisten zum System beitragenden Ordnungen werden getrackt und
% dazugehörend magnitude, rpm und time bestimmt.
[mag, rpm, ~] = ordertrack(vibData, fs, [rpmVec rpmVec], orderlist, [1 2]);

mag = mag';
orderlist = orderlist';
rpm = rpm';

% Frequenzband der jeweiligen orders
freqLine = rpm/60.*orderlist(:)';

mag_norm = mag(68:end,:) ./ max(mag(68:end,:),[],1);

% 2. Subplot. Visualisiert das Frequenzspektrum
subplot(1,2,2)
plot(freqLine(68:end,:), mag_norm);

end

