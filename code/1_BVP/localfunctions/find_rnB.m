function rnB = find_rnB(sig,fs,vnA,window)

% sig       time signal @ fs
% fs        sampling freq [Hz]
% nA        peak annotations [samples]
% window    backward shift for search window nB

% rnB       local min in sig(nA-window:nA) with edge correction


%% note: edge cases might need to be cleaned in RDECO!


%% sanity check

try
    assert(window<5000/fs) % prob. 150ms 
catch
    disp('Most likely search window too large, reconsider param')
end


%% Apply basic search algo for each nA

rnB = zeros(size(vnA));

for i = 1:length(vnA)

    if vnA(i)-window+1 <= 1
        start=1;
    else
        start=vnA(i)-window+1;
    end
    
    local = sig(uint16(start:vnA(i)));
    [~,ref_idx] = min(local);
    
    rnB(i) = vnA(i) - window + ref_idx;
    
end


end