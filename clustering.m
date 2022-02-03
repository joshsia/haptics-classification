% Section D

clear all; close all; clc;

load('data-matrices/F0_PVT.mat')
load('data-matrices/F0_Electrodes.mat')

all_objects{1} = acrylic_pvt;
all_objects{2} = black_foam_pvt;
all_objects{3} = car_sponge_pvt;
all_objects{4} = flour_sack_pvt;
all_objects{5} = kitchen_sponge_pvt;
all_objects{6} = steel_vase_pvt;

all_electrodes = {acrylic_electrodes,black_foam_electrodes,...
    car_sponge_electrodes,flour_sack_electrodes,...
    kitchen_sponge_electrodes,steel_vase_electrodes};

input_pvt = zeros(60,3);
input_elec = zeros(60,19);
idx = 1:10:60;
for i = 1:6
input_pvt(idx(i):idx(i)+9,1) = all_objects{i}.pressure;
input_pvt(idx(i):idx(i)+9,2) = all_objects{i}.vibrations;
input_pvt(idx(i):idx(i)+9,3) = all_objects{i}.temp;
input_elec(idx(i):idx(i)+9,:) = all_electrodes{i}';
end
% input_pvt = zscore(input_pvt);
% input_elec = zscore(input_elec);

y = [ones(1,10) 2*ones(1,10) 3*ones(1,10) 4*ones(1,10) 5*ones(1,10) 6*ones(1,10)];
y_label = cell(1,60);
idx = 1:10:60;
for i = 1:6
    for j = idx(i) : idx(i)+9
        y_label{j} = num2str(i);
    end
end

plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#A2142F"];
plot_colours2 = [0,114,189; 217,83,25; 237,177,32; 126,47,142; 119,172,48; 162,20,47]/255;


% Clustering PVT

rng(223);

metric = "seuclidean";
% metric = "mahalanobis";
tree = linkage(input_pvt,'average',metric);
figure();
cutoff = median([tree(end-5,3) tree(end-4,3)]);
P = 0;
H = dendrogram(tree,P,'ColorThreshold',cutoff,'Labels',y_label);
% H = dendrogram(tree,P,'ColorThreshold',cutoff);

lineColours = cell2mat(get(H,'Color'));
colourList = unique(lineColours, 'rows');

idxs = zeros(1,size(colourList,1)-1);
for colour = 2:size(colourList,1)
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    idxs(colour-1) = sum(idx);
    lineColours(idx, :) = repmat(plot_colours2(colour-1,:),sum(idx),1);
end

for line = 1:size(H,1)
    set(H(line),'Color',lineColours(line,:));
    H(line).LineWidth = 1.0;
end

ax = gca;
ax.FontSize = 16; 
ax.XAxis.FontSize = 12;
title(sprintf('Dendogram of the PVT data based on the %s metric',metric),'FontSize',13)
xlabel('True class','FontSize',16)
ylabel('Distance','FontSize',16)
% saveas(gcf,sprintf('Images/Section D/pvt_cluster_dendogram_%s.eps',metric),'epsc');

T = cluster(tree,'maxclust',6);
% crosstab(T,y)

my_classes = unique(T);
n_each = zeros(2,length(my_classes));
for i = 1:length(my_classes)
    n_each(1,i) = sum(T==i)-1;
    n_each(2,i) = i;
end

for i = 1:length(my_classes)
    for j = i+1:length(my_classes)
        if(n_each(1,j)==idxs(i))
            temp = n_each(:,j);
            n_each(:,j) = n_each(:,i);
            n_each(:,i) = temp;
        end
    end
end

T_colour = zeros(size(T));
for i = 1:length(T)
    T_colour(i) = find(n_each(2,:)==T(i));
end

fig = figure('Position', [100, 0, 1000, 900]);
subplot(1,5,[1 4])
for i = 1:length(T)
    scatter3(input_pvt(i,1),...
        input_pvt(i,2),...
        input_pvt(i,3),75,...
        'o','filled','MarkerFaceColor',plot_colours(T_colour(i)))
    hold on;
end
ax = gca;
ax.FontSize = 16; 
xlabel('Pressure','FontSize',16)
ylabel('Vibrations','FontSize',16)
zlabel('Temperature','FontSize',16)
title(sprintf('Clusters found by hierarchical clustering of the PVT data (%s)',metric))

hSub = subplot(1,5,5);
hold on
for i = 1:6
    plot([1 nan],[1 nan],'Color',plot_colours(i),'LineWidth',1.0)
end
set(hSub, 'Visible', 'off');
legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4',...
    'Cluster 5','Cluster 6'},...
    'FontSize',16,'Location','east')
% saveas(gcf,sprintf('results/clustering/pvt_cluster_3dplot_%s.eps',metric),'epsc');

%%
% Classification of Electrodes data
load('data-matrices/Electrodes_PCA.mat')

data = score_3D;
y = repmat([1:1:6],10,1);
y = y(:);
X = [data y];

frac = 1;
nTrees = 1:1:100;

% Loop to find the optimal number of trees to be used
% seeds = [223 123 456 2021 2000 178 345 999 1 9];
% all_acc = zeros(length(seeds),length(nTrees));
% for my_seed = 1:length(seeds)
%     rng(seeds(my_seed));
%     
%     cv = cvpartition(size(X,1),'HoldOut',0.4);
%     idx = cv.test;
%     train_X = X(~idx,:);
%     test_X = X(idx,:);
% 
%     for i = 1:length(nTrees)
%         model = TreeBagger(nTrees(i),train_X(:,1:3),train_X(:,4),'Method','classification',...
%             'InBagFraction',frac);
% 
%         temp = predict(model,test_X(:,1:3));
% 
%         pred = zeros(length(temp),1);
%         for j = 1:length(temp)
%             pred(j) = str2double(temp{j});
%         end
% 
%         all_acc(my_seed,i) = 100*sum(test_X(:,4)==pred)/length(pred);
%     end
% end
% 
% mean_acc = mean(all_acc,1);
% opt_nTrees = find(mean_acc==max(mean_acc),1);

% figure()
% plot(nTrees, mean_acc);
% ylabel('Test Accuracy','FontSize',16);
% xlabel('Number of trees','FontSize',16);
% title('Electrode: Test Accuracy as function of Number of Trees','FontSize',14);
% saveas(gcf,'results/clustering/elec_mean_acc.eps','epsc');

% Optimal number of trees was found to be 18
rng(178);
opt_nTrees = 18;

cv = cvpartition(size(X,1),'HoldOut',0.4);
idx = cv.test;
train_X = X(~idx,:);
test_X = X(idx,:);

model = TreeBagger(opt_nTrees,train_X(:,1:3),train_X(:,4),'Method','classification',...
 'InBagFraction',frac,'OOBPrediction','on');
% save('elecClassifyModel.mat','model')
% model = load('elecClassifyModel.mat');
err = oobError(model);

figure()
plot(nTrees(1:18), err);
xlabel('Number of Trees','FontSize',16)
ylabel('Out-of-bag error','FontSize',16)
title('Electrode: Out-of-bag error against Number of Trees','FontSize',14)
% saveas(gcf,'results/clustering/elec_oob.eps','epsc');

temp = predict(model,test_X(:,1:3));

pred = zeros(length(temp),1);
for j = 1:length(temp)
    pred(j) = str2double(temp{j});
end
    
view(model.Trees{1},'mode','graph')
view(model.Trees{2},'mode','graph')

C = confusionmat(test_X(:,4),pred);
labels = ["Acrylic", "Black foam", "Car sponge","Flour sack",...
    "Kitchen sponge","Steel vase"];

figure();
cc = confusionchart(C,labels,'FontSize',13);
cc.title('Electrode: Confusion Chart');
% saveas(gcf,'results/clustering/elec_confusion_chart.eps','epsc');


