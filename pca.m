% Section B

clear all; close all; clc;

trials = 10;
n_objects = 6;

load('data-matrices/F0_PVT.mat')
load('data-matrices/F0_Electrodes.mat')

acrylic = [acrylic_pvt.vibrations;acrylic_pvt.pressure;acrylic_pvt.temp];
black_foam = [black_foam_pvt.vibrations;black_foam_pvt.pressure;black_foam_pvt.temp];
car_sponge = [car_sponge_pvt.vibrations;car_sponge_pvt.pressure;car_sponge_pvt.temp];
flour_sack = [flour_sack_pvt.vibrations;flour_sack_pvt.pressure;flour_sack_pvt.temp];
kitchen_sponge = [kitchen_sponge_pvt.vibrations;kitchen_sponge_pvt.pressure;kitchen_sponge_pvt.temp];
steel_vase = [steel_vase_pvt.vibrations;steel_vase_pvt.pressure;steel_vase_pvt.temp];

pvt_data = [acrylic black_foam car_sponge flour_sack kitchen_sponge steel_vase]';
electrode_data = [acrylic_electrodes black_foam_electrodes car_sponge_electrodes ...
    flour_sack_electrodes kitchen_sponge_electrodes steel_vase_electrodes]';

plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#000000"];

% PCA on PVT data
Z = zscore(pvt_data);
[coefs,score,latent] = pca(Z);

figure();
vbls = {'Vibrations','Pressure','Temperature'};
h1 = biplot(coefs(:,1:3),'Scores',score(:,1:3),'VarLabels',vbls,...
    'MarkerSize',12);
ax = gca;
ax.FontSize = 16; 
title('Standardised PVT data with 3 Principal Components','FontSize',14)
xlabel('Component 1','FontSize',16)
ylabel('Component 2','FontSize',16)
zlabel('Component 3','FontSize',16)
% saveas(gcf, 'results/pca/pvt_pca_3d.eps','epsc')


for k = 1:6
    h1(k).Color = 'b';
    h1(k).LineWidth = 1.5;
end

for k = 7:9
    h1(k).FontWeight = 'bold';
    h1(k).FontSize = 14;
end

obj = [10 20 30 40 50 60];
for i = 1:length(obj)
    for k = obj(i) : obj(i)+10
        h1(k).MarkerEdgeColor = plot_colours(i);
    end
end

% 2D reduction
coefs_2D = coefs(:,1:2);
score_2D = score(:,1:2);

figure();
vbls = {'Vibrations','Pressure','Temperature'};
h2 = biplot(coefs_2D,'Scores',score_2D,'VarLabels',vbls,'MarkerSize',12);
ax = gca;
ax.FontSize = 16; 
title('Standardised PVT data with 2 Principal Components','FontSize',14)
xlabel('Component 1','FontSize',16)
ylabel('Component 2','FontSize',16)
% saveas(gcf, 'results/pca/pvt_pca_2d.eps','epsc')


for k = 1:6
    h2(k).Color = 'b';
    h2(k).LineWidth = 1.5;
end

for k = 7:9
    h2(k).FontWeight = 'bold';
    h2(k).FontSize = 14;
end

obj = [10 20 30 40 50 60];
for i = 1:length(obj)
    for k = obj(i) : obj(i)+10
        h2(k).MarkerEdgeColor = plot_colours(i);
    end
end

temp_axis = [-0.52593 0.84501];
vibrations_axis = [-0.59439 -0.44649];
pressure_axis = [0.60836 0.29428];

% Pressure smallest magnitude
temp_mag = sum(temp_axis.^2);
vib_mag = sum(vibrations_axis.^2);
press_mag = sum(pressure_axis.^2);

PC = Z*coefs;

features = length(latent);
figure();
for i = 1:features
    % subplot(features,1,i);
    for j = 1:length(PC)
        hold on
        col = 1 + floor((j-1)/10);
        plot(PC(j,i),i,'.','MarkerSize',15,'Color',plot_colours(col));
    end
   
end
xlim([-2.5 2.5])
ylim([0 4])
yticks([1 2 3])
set (gca,'ydir','reverse');
ax = gca;
ax.FontSize = 16; 
title('Distribution of the Standardised PVT data along 3 Principal Components','FontSize',14)
xlabel('Score','FontSize',16)
ylabel('Component number','FontSize',16)
% saveas(gcf, 'results/pca/pvt_pca_1d.eps','epsc')


% PCA function on electrodes
Z = zscore(electrode_data);
[coefs,score,latent] = pca(Z);

% Scree plot (variance)
figure();
plot(latent,'LineWidth',1.5)
ax = gca;
ax.FontSize = 16; 
title('Scree plot:  Variance against Principal Component number for electrode data','FontSize',14)
xlabel('Component number','FontSize',16)
ylabel('Variance','FontSize',16)
ylim([-0.5 13])
grid on
grid minor
% saveas(gcf, 'results/pca/var_scree_plot.eps','epsc')

% Scree plot (log variance)
% figure();
% plot(log(latent),'LineWidth',1.5)
% ax = gca;
% ax.FontSize = 16; 
% title('Electrode: Log Variance against Principal Components','FontSize',14)
% xlabel('Component number','FontSize',16)
% ylabel('Log variance','FontSize',16)
% grid on
% grid minor
% % saveas(gcf, 'results/pca/logvar_scree_plot.eps','epsc')

% Reduce to 3D
coefs_3D = coefs(:,1:3);
score_3D = score(:,1:3);
latent_3D = latent(1:3);

figure();
h3 = biplot(coefs_3D,'Scores',score_3D,'MarkerSize',12);
ax = gca;
ax.FontSize = 16; 
title('Standardised Electrode data with 3 Principal Components','FontSize',14)
xlabel('Component 1','FontSize',16)
ylabel('Component 2','FontSize',16)
zlabel('Component 3','FontSize',16)
% saveas(gcf, 'results/pca/elec_pca_3d.eps','epsc')

for k = 1:19
    h3(k).Color = 'b';
    h3(k).LineWidth = 1.0;
end

for k = 39:length(h3)
    h3(k).MarkerEdgeColor = 'k';
end


PC = Z*coefs;

features = length(latent);
figure();
for i = 1:features
    for j = 1:length(PC)
        hold on
        col = 1 + floor((j-1)/10);
        plot(PC(j,i),i,'.','MarkerSize',15,'Color','k');
    end
    
end
ylim([0 20])
yticks([1 2 3 4 5 6 7 8 9 10 ...
    11 12 13 14 15 16 17 18 19])
set (gca,'ydir','reverse');
ax = gca;
ax.FontSize = 16; 
title({'Distribution of the Standardised Electrode data','along 19 Principal Components'},'FontSize',14)
xlabel('Score','FontSize',16)
ylabel('Component number','FontSize',16)
% saveas(gcf, 'results/pca/elec_pca_1d.eps','epsc')
% save('data-matrices/Electrodes_PCA.mat','Z','coefs_3D','score_3D','latent_3D','PC');
