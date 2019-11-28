clc
clear all

%% Divide indentation sessions
k1=1;
k1_old=0;
n=3;
for probe = 1 : 4
    if probe == 1 
        name = 'A';
        elseif probe == 2
            name = 'B';
        elseif probe == 3
            name = 'C';
        elseif probe == 4
            name = 'D';
            n=2;
    end
    
for i =0 : n-1
if i==0 
    s = strcat('./dataset/test_16channels_Probe', name,'.xlsx');
else
    ind = int2str(i);
    s = strcat('./dataset/test_16channels_Probe', name,'_', ind, '.xlsx');
end

Data = readtable(s);
Data = table2array(Data);
Data(:,5)= probe;

ne0 = find(diff(Data(:,1))> 300);

section{k1} = Data (1:ne0(1),:);
for j = k1+1:length(ne0)+k1_old
    
    section{j} = Data(ne0(j-1-k1_old)+1:ne0(j-k1_old),:);
end
%section{k1+1} = Data(ne0(k1-k1_old):length(Data),:);
k1= k1 + length(ne0);
k1_old = k1-1;
end
end
section(1) = [];
section(132) = [];
save('dataset','Data','section')

