function windowsize = get_valid_windowsize(windowlength,fs,up)

% windowlength      duration of given window [sec]
% fs                sampling freq for used signal [Hz]
% up                round up if True [logical]

% windowsize        size of valid window [samples]: even integer


%% Get nearest even integer

size_raw = windowlength*fs;
nearest_int = round(size_raw);

if mod(nearest_int,2)==1
    if up
        size_valid = nearest_int + 1;
    else
        size_valid = nearest_int - 1;
    end
else
    size_valid = nearest_int;
end

windowsize = size_valid;

end