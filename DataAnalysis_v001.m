%% RUN AFTER DataElaboration_v002.m 

%% Spike Count

% for each indentation session in "section" counts the spike for each channel. 
% section{indentation}(:,2)==channel returns a vector of 1 if the channel
% selected spikes 0 otherwise by summing these elements we obtain the spike
% count

for indentation = 1:size(section,2)
    for channel = 1:16
        spike_count(indentation,channel) = sum(section{indentation}(:,2)==channel);
    end
    if spike_count(indentation,i) == spike_count(indentation,2) & spike_count(indentation,2) == spike_count(indentation,3) & spike_count(indentation,3) == spike_count(indentation,4)
        spike_count(indentation,:) = []; 
        section(indentation)=[];
    end
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
distance = zeros (size(section,2),4);

% for each session computes the VP distance on each channel and sums them
% to obtain a cumulative one
for indentation = 1 :size(section,2)
    data = zeros(21,16); 
    for channel = 1:16
        index = section{indentation}(:,2)==channel;
        data(1:sum(index),channel) = section{indentation}(index,1);
        data(1:sum(index),channel) = data(1:sum(index),channel) - data(1,channel);
        distance(indentation,1) = distance(indentation,1) + victor_purpura_distance_f(data(:,channel),data_A(:,channel),0.01);
        distance(indentation,2) = distance(indentation,2) + victor_purpura_distance_f(data(:,channel),data_B(:,channel),0.01);
        distance(indentation,3) = distance(indentation,3) + victor_purpura_distance_f(data(:,channel),data_C(:,channel),0.01);
        distance(indentation,4) = distance(indentation,4) + victor_purpura_distance_f(data(:,channel),data_D(:,channel),0.01);
    end
end

% distance is a vector that has probe type as columns and indentation no.
% as row

%% Plot results

indNo = 1:size(distance,1);
subplot(2,2,1)
hold on
patch([0 66 66 0], [65 65 0 0], 'r')
patch([66 114 114 66], [65 65 0 0], 'g')
patch([114 182 182 114], [65 65 0 0], 'c')
patch([182 202 202 182], [65 65 0 0], 'm')
alpha(0.3)
ylim([0 65]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,1), 'b.')
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
alpha(0.3)
ylim([0 65]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,2), 'b.')
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
alpha(0.3)
ylim([0 65]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,3), 'b.')
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
alpha(0.3)
ylim([0 65]);
xlim([0 size(distance,1)]);
plot( indNo,distance(:,4), 'b.')
ylabel('Victor Purpura distance')
xlabel('Indentation session')
title('Distance from probe D')
hold off
