function rnM = find_rnM(sig,vnA,vnB)

% sig       time signal @ fs
% vnA       validated/clean apex point annotations [samples]
% vnB       validated/clean basepoint annotations [samples]

% rnM       midpoint annotations, 50% amplitude sig(nB:nA)


%% Remarks

% % Assumptions
% every nB matches a nA, in order of nB before nA
% edge case errors and artefacts should have been removed/solved

% % Code reference
% algo from Jesus Lazaro, ComputePPGPeakPoints.m line 301


%% sanity checks

try
    assert(all(size(vnA)==size(vnB)));
catch
    disp('WARNING - annotations do not match size, check edge cases')
end

try
    assert(isempty(find(vnA-vnB<0,1)));
catch
    disp('WARNING - vnA < vnB for 1+ annotations, check input args')
end


%% Apply algo (ref003, Lazaro) for each detected pulse

rnM = zeros(size(vnA));

for i = 1:length(vnA)

    % Exact mean amplitude of signal between nB and nA:
    mean_amp = 0.5*(sig(vnB(i))+sig(vnA(i)));

    % Find sample closest to exact mean amplitude in the range [nB..nA]
    searchlist = vnB(i):vnA(i);
    amplitudes = sig(searchlist);
    [~,aux_pos]=min(abs(amplitudes-mean_amp)); 
    rnM(i)=vnB(i)+aux_pos-1;

end


end