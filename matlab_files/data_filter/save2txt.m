%%
load('t.mat');
fid = fopen('t.txt','w');
row = length(t);
for i=1:row
    fprintf(fid,'%g%s',t(i),' ');
    if i~=row
        fprintf(fid,'\n');
    end
end
fclose(fid);

load('q_filter.mat');
fid = fopen('q_filter.txt','w');
[m,n] = size(q_filter);
for i=1:m
    for j=1:n
        fprintf(fid,'%g%s',q_filter(i,j),' ');
    end
    if i~=m
        fprintf(fid,'\n');
    end
end
fclose(fid);

load('qDot_filter.mat');
fid = fopen('qDot_filter.txt','w');
[m,n] = size(qDot_filter);
for i=1:m
    for j=1:n
        fprintf(fid,'%g%s',qDot_filter(i,j),' ');
    end
    if i~=m
        fprintf(fid,'\n');
    end
end
fclose(fid);

load('qDDot_filter.mat');
fid = fopen('qDDot_filter.txt','w');
[m,n] = size(qDDot_filter);
for i=1:m
    for j=1:n
        fprintf(fid,'%g%s',qDDot_filter(i,j),' ');
    end
    if i~=m
        fprintf(fid,'\n');
    end
end
fclose(fid);

load('tau_filter.mat');
fid = fopen('tau_filter.txt','w');
[m,n] = size(tau_filter);
for i=1:m
    for j=1:n
        fprintf(fid,'%g%s',tau_filter(i,j),' ');
    end
    if i~=m
        fprintf(fid,'\n');
    end
end
fclose(fid);
