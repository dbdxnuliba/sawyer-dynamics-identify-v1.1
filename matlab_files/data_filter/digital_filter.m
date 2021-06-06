function data_filter = digital_filter(n, fs, fc, data)

[b,a] = butter(n,fc/(fs/2));
data_filter=filtfilt(b, a, data);

end