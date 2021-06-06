function data_diff = digital_differentiator(fs, data)

Ts=1/fs;
n=length(data);
data_diff=zeros(n,1);

for k=2:(n-1)
    data_diff(k,1)=(data(k+1,1)-data(k-1,1))/(2*Ts);
end

data_diff(1,1)=data_diff(2,1);
data_diff(n,1)=data_diff(n-1,1);

end