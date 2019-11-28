close all
clear all

load('dataset')

%% Spike Count

% for each indentation session in "section" counts the spike for each channel. 
% section{indentation}(:,2)==channel returns a vector of 1 if the channel
% selected spikes 0 otherwise by summing these elements we obtain the spike
% count

for indentation = 1:size(section,2)
    for channel = 1:16
        spike_count(indentation,channel) = sum(section{indentation}(:,2)==channel);
    end
    spike_no(indentation) = size(section{indentation},1);
end

%% Victor Purpura distance

      data_A = zeros(21,16); 
      data_B = zeros(21,16); 
      data_C = zeros(21,16); 
      data_D = zeros(21,16);
      
 % Reorganize datas so that every indentation session has channels as columns
 % spike count as rows and spike times as elements. Takes 4 known sessions
 % A, B, C, D from which measure the Victor Purpura distance 
 
 
% reorganize the 4 known sessions
for channel = 1:16 
      index_A = section{1}(:,2)==channel;
      data_A(1:sum(index_A),channel) = section{1}(index_A,1);
      data_A(1:sum(index_A),channel) = data_A(1:sum(index_A),channel) - data_A(1,channel);
      index_B = section{67}(:,2)==channel;
      data_B(1:sum(index_B),channel) = section{67}(index_B,1);
      data_B(1:sum(index_B),channel) = data_B(1:sum(index_B),channel) - data_B(1,channel);
      index_C = section{114}(:,2)==channel;
      data_C(1:sum(index_C),channel) = section{114}(index_C,1);
      data_C(1:sum(index_C),channel) = data_C(1:sum(index_C),channel) - data_C(1,channel);
      index_D = section{186}(:,2)==channel;
      data_D(1:sum(index_D),channel) = section{186}(index_D,1);
      data_D(1:sum(index_D),channel) = data_D(1:sum(index_D),channel) - data_D(1,channel);
end
distance = zeros (size(section,2),5);
null =  zeros(21,16); 

% for each session computes the VP distance on each channel and sums them
% to obtain a cumulative one
for indentation = 1 :size(section,2)
    data = zeros(21,16);
    for channel = 1:16
        index = section{indentation}(:,2)==channel;
        data(1:sum(index),channel) = section{indentation}(index,1) + rand(sum(index),1)*0;
        data(1:sum(index),channel) = data(1:sum(index),channel) - data(1,channel); 
        distance(indentation,1) = distance(indentation,1) + victor_purpura_distance_f(data(:,channel),data_A(:,channel),0.001);
        distance(indentation,2) = distance(indentation,2) + victor_purpura_distance_f(data(:,channel),data_B(:,channel),0.001);
        distance(indentation,3) = distance(indentation,3) + victor_purpura_distance_f(data(:,channel),data_C(:,channel),0.001);
        distance(indentation,4) = distance(indentation,4) + victor_purpura_distance_f(data(:,channel),data_D(:,channel),0.001);
        distance(indentation,5) = distance(indentation,5) + victor_purpura_distance_f(data(:,channel),null(:,channel),0.001);
    end
    if section{indentation}(1,5) == 1
        label{indentation,1} = 'A';
    end
    if section{indentation}(1,5) == 2
        label{indentation,1} = 'B';
    end
    if section{indentation}(1,5) == 3
        label{indentation,1} = 'C';
    end
    if section{indentation}(1,5) == 4
        label{indentation,1} = 'D';
    end
end

% distance is a vector that has probe type as columns and indentation no.
% as row

DistanceCountA = zeros(size(spike_count,1),2);
DistanceCountA(:,1) = spike_count(:,9);
DistanceCountA(:,2) = distance(:,1);
DistanceCountB = zeros(size(spike_count,1),2);
DistanceCountB(:,1) = spike_count(:,9);
DistanceCountB(:,2) = distance(:,2);
DistanceCountC = zeros(size(spike_count,1),2);
DistanceCountC(:,1) = spike_count(:,9);
DistanceCountC(:,2) = distance(:,3);
DistanceCountD = zeros(size(spike_count,1),2);
DistanceCountD(:,1) = spike_count(:,9);
DistanceCountD(:,2) = distance(:,4);
DistanceCountNull = zeros(size(spike_count,1),2);
DistanceCountNull(:,1) = spike_count(:,9);
DistanceCountNull(:,2) = distance(:,5);


%% Clustering KKN

centroids = zeros(3,2);
centroids(1,:) = DistanceCountNull(1,:);
centroids(2,:) = DistanceCountNull(67,:);
centroids(3,:) = DistanceCountNull(114,:);
centroids(4,:) = DistanceCountNull(186,:);
percentage = round(0.7*size(DistanceCountNull,1));
s = RandStream('mt19937ar','Seed',0);
indices(:,1) = randperm(s,size(DistanceCountNull,1));
indices_train = indices(1:percentage);
indices_test = indices(percentage+1:end);
train_set = DistanceCountNull(indices_train,:);
test_set = DistanceCountNull(indices_test,:);

lda = fitcknn(train_set,label(indices_train));
ldaClass = resubPredict(lda);
ldaResubErr = resubLoss(lda);
test_pred = predict(lda, test_set);
bad = ~ strcmp(ldaClass,label(indices_train));
x1range = min(DistanceCountNull(:,1))-2:.01:max(DistanceCountNull(:,1))+2;
x2range = min(DistanceCountNull(:,2))-2:.01:max(DistanceCountNull(:,2))+2;
[xx1, xx2] = meshgrid(x1range,x2range);
XGrid = [xx1(:) xx2(:)];
predictedspecies = predict(lda,XGrid);

%% Plot results

indNo = 1:size(distance,1);
th = ones(size(distance,1))*2;

fig2=figure();
subplot(2,2,1)
hold on
patch([0 66 66 0], [65 65 0 0], 'r')
patch([66 114 114 66], [65 65 0 0], 'g')
patch([114 182 182 114], [65 65 0 0], 'c')
patch([182 202 202 182], [65 65 0 0], 'm')
alpha(0.5)
ylim([0 8]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,1), 'b.', 'MarkerSize',8)
plot( indNo,th, 'w--', 'LineWidth', 2)
ylabel('Victor Purpura distance')
xlabel('Indentation session')
title('Distance from probe A')
hold off

subplot(2,2,2)
hold on
patch([0 66 66 0], [65 65 0 0], 'r')
patch([66 114 114 66], [65 65 0 0], 'g')
patch([114 182 182 114], [65 65 0 0], 'c')
patch([182 202 202 182], [65 65 0 0], 'm')
alpha(0.5)
ylim([0 8]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,2), 'b.', 'MarkerSize',8)
plot( indNo,th, 'w--', 'LineWidth', 2)
ylabel('Victor Purpura distance')
xlabel('Indentation session')
title('Distance from probe B')
hold off

subplot(2,2,3)
hold on
patch([0 66 66 0], [65 65 0 0], 'r')
patch([66 114 114 66], [65 65 0 0], 'g')
patch([114 182 182 114], [65 65 0 0], 'c')
patch([182 202 202 182], [65 65 0 0], 'm')
alpha(0.5)
ylim([0 8]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,3), 'b.', 'MarkerSize',8)
plot( indNo,th, 'w--', 'LineWidth', 2)
ylabel('Victor Purpura distance')
xlabel('Indentation session')
title('Distance from probe C')
hold off

subplot(2,2,4)
hold on
patch([0 66 66 0], [65 65 0 0], 'r')
patch([66 114 114 66], [65 65 0 0], 'g')
patch([114 182 182 114], [65 65 0 0], 'c')
patch([182 202 202 182], [65 65 0 0], 'm')
alpha(0.5)
ylim([0 8]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,4), 'b.', 'MarkerSize',8)
plot( indNo,th, 'w--', 'LineWidth', 2)
ylabel('Victor Purpura distance')
xlabel('Indentation session')
title('Distance from probe D')
hold off

fig1 = figure();
alpha(.3);
f1 = gscatter(xx1(:), xx2(:), predictedspecies,'rgcm');
ylim([min(DistanceCountNull(:,2))-1 max(DistanceCountNull(:,2))+1]);
xlim([min(DistanceCountNull(:,1))-1 max(DistanceCountNull(:,1))+1]);
hold on
plot(centroids(:,1),centroids(:,2),'k*')
gscatter(DistanceCountNull(:,1),DistanceCountNull(:,2), label, 'brwc', 'osdp')
plot(DistanceCountNull(bad,1),DistanceCountNull(bad,2), 'rx')
legend('Cluster A','Cluster B','Cluster C','Cluster D')
title('Clustering')
xlabel('Spike count on channnel 9')
ylabel('Victor purpura distance')


figure;
subplot(1,2,1)
ldaResubCM = confusionchart(label(indices_train),ldaClass);
title('Performance on Train-set')
subplot(1,2,2)
confusionchart(label(indices_test),test_pred);
title('Performance on Test-set')

