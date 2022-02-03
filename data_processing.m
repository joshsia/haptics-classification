% Section A - 1

% 'F0Electrodes','F1Electrodes', - Electrode Impedance
% 'F0pac','F1pac' - High Frequency Fluid Vibrations
% 'F0pdc','F1pdc' - Low Frequency Fluid Pressure
% 'F0tac','F1tac' - Core Temperature Change
% 'F0tdc','F1tdc' - Core Temperature
% 'JEff' - Robot arm joint effort (load)
% 'JPos' � Robot arm joint positions
% 'JVel' - Robot arm joint velocity

close all; clear;

% Loading datasets
trials = 10;
n_objects = 6;
objects = cell(trials,n_objects);

for i = 1:trials
    a_temp = "data/acrylic_211_0" + num2str(i) + "_HOLD.mat";
    b_temp = "data/black_foam_110_0" + num2str(i) + "_HOLD.mat";
    c_temp = "data/car_sponge_101_0" + num2str(i) + "_HOLD.mat";
    f_temp = "data/flour_sack_410_0" + num2str(i) + "_HOLD.mat";
    k_temp = "data/kitchen_sponge_114_0" + num2str(i) + "_HOLD.mat";
    s_temp ="data/steel_vase_702_0" + num2str(i) + "_HOLD.mat";
    if(i==10)
        a_temp = "data/acrylic_211_10_HOLD.mat";
        b_temp = "data/black_foam_110_10_HOLD.mat";
        c_temp = "data/car_sponge_101_10_HOLD.mat";
        f_temp = "data/flour_sack_410_10_HOLD.mat";
        k_temp = "data/kitchen_sponge_114_10_HOLD.mat";
        s_temp = "data/steel_vase_702_10_HOLD.mat";
    end
    objects{i,1} = load(a_temp);
    objects{i,2} = load(b_temp);
    objects{i,3} = load(c_temp);
    objects{i,4} = load(k_temp);
    objects{i,5} = load(f_temp);
    objects{i,6} = load(s_temp);
end

% Finding mean values
% i==2 and j==8 is of length 1x1001. Take only 1:1000

mean_objects = cell(1,n_objects);
for i = 1:n_objects
    temp_pac = zeros(1,1000);
    temp_pdf = zeros(1,1000);
    temp_tdc = zeros(1,1000);
    temp_electrodes = zeros(19,1000);
    for j = 1:trials
        temp_pac = temp_pac + (1/trials).*objects{j,i}.F0pac(2,1:1000);
        temp_pdf = temp_pdf + (1/trials).*objects{j,i}.F0pdc(1:1000);
        temp_tdc = temp_tdc + (1/trials).*objects{j,i}.F0tdc(1:1000);
        temp_electrodes = temp_electrodes + (1/trials).*objects{j,i}.F0Electrodes(:,1:1000);
    end
    mean_objects{i} = struct('vibrations',temp_pac,...
        'pressure',temp_pdf,...
        'temperature',temp_tdc,...
        'electrodes',temp_electrodes);
end

% Finding time step that allows differentiation between different objects

% Method: Given a score based on (i) Inter-object variance at time step,
% (ii) Smallest distance between two objects at time step,
% (iii) Intra-object variance at time step

standard_pac = zeros(10,6,1000);
standard_pdf = zeros(10,6,1000);
standard_tdc = zeros(10,6,1000);
for i = 1:trials
    for j = 1:n_objects
        standard_pac(i,j,:) = objects{i,j}.F0pac(2,1:1000);
        standard_pdf(i,j,:) = objects{i,j}.F0pdc(1:1000);
        standard_tdc(i,j,:) = objects{i,j}.F0tdc(1:1000);
    end
end

standard_pac = zscore(standard_pac,0,[1 2]);
standard_pdf = zscore(standard_pdf,0,[1 2]);
standard_tdc = zscore(standard_tdc,0,[1 2]);

all_sd_pac = zeros(trials,1000);
all_sd_pdf = zeros(trials,1000);
all_sd_tdc = zeros(trials,1000);
min_diff_pac = zeros(trials,1000);
min_diff_pdf = zeros(trials,1000);
min_diff_tdc = zeros(trials,1000);
for i = 1:trials
    temp_pac = zeros(n_objects,1000);
    temp_pdf = zeros(n_objects,1000);
    temp_tdc = zeros(n_objects,1000);
    for j = 1:n_objects
        temp_pac(j,:) = standard_pac(i,j,:);
        temp_pdf(j,:) = standard_pdf(i,j,:);
        temp_tdc(j,:) = standard_tdc(i,j,:);
    end
    all_sd_pac(i,:) = std(temp_pac);
    all_sd_pdf(i,:) = std(temp_pdf);
    all_sd_tdc(i,:) = std(temp_tdc);
    
    for k = 1:1000
        min_diff = 1/eps;
        min_pdf = 1/eps;
        min_tdc = 1/eps;
        for l = 1:n_objects
            for m = (l+1):n_objects
                temp_diff = (temp_pac(l,k) - temp_pac(m,k)).^2;
                if(temp_diff < min_diff)
                    min_diff = temp_diff;
                end
                temp_diff = (temp_pdf(l,k) - temp_pdf(m,k)).^2;
                if(temp_diff < min_pdf)
                    min_pdf = temp_diff;
                end
                temp_diff = (temp_tdc(l,k) - temp_tdc(m,k)).^2;
                if(temp_diff < min_tdc)
                    min_tdc = temp_diff;
                end
            end
        end
        min_diff_pac(i,k) = min_diff;
        min_diff_pdf(i,k) = min_pdf;
        min_diff_tdc(i,k) = min_tdc;
    end
end

mean_sd_pac = mean(all_sd_pac,1);
mean_sd_pdf = mean(all_sd_pdf,1);
mean_sd_tdc = mean(all_sd_tdc,1);
mean_min_diff_pac = mean(min_diff_pac,1);
mean_min_diff_pdf = mean(min_diff_pdf,1);
mean_min_diff_tdc = mean(min_diff_tdc,1);

intra_pac = zeros(n_objects,1000);
intra_pdf = zeros(n_objects,1000);
intra_tdc = zeros(n_objects,1000);
for i = 1:n_objects
    temp_pac = zeros(trials,1000);
    temp_pdf = zeros(trials,1000);
    temp_tdc = zeros(trials,1000);
    for j = 1:trials
        temp_pac(j,:) = standard_pac(j,i,:);
        temp_pdf(j,:) = standard_pdf(j,i,:);
        temp_tdc(j,:) = standard_tdc(j,i,:);
    end
    intra_pac(i,:) = std(temp_pac);
    intra_pdf(i,:) = std(temp_pdf);
    intra_tdc(i,:) = std(temp_tdc);
end

total_intra_pac = sum(intra_pac,1);
total_intra_pdf = sum(intra_pdf,1);
total_intra_tdc = sum(intra_tdc,1);

% We want min_diff and sd to be as big as possible
% We want intra-object variance to be small
% Give each time step a score

% Gives time step = 38
weights = [0.35 0.45 0.2];
score = weights(1)*(mean_sd_pac + mean_sd_pdf + mean_sd_tdc) + ...
    weights(2)*(mean_min_diff_pac + mean_min_diff_pdf + mean_min_diff_tdc) + ...
    + weights(3)*(- total_intra_pac - total_intra_pdf - total_intra_tdc);

% % Gives time step = 695
% score = (mean_min_diff_pac + mean_min_diff_pdf + mean_min_diff_tdc);

time_step = find(score==max(score));

% figure();
% plot(score)

% figure();
% subplot(311)
% plot(all_sd_pac.')
% ax = gca;
% ax.FontSize = 16; 
% ylabel('Vibrations','FontSize',16)
% xlabel('Time','FontSize',16)
% subplot(312)
% plot(all_sd_pdf.')
% ax = gca;
% ax.FontSize = 16; 
% ylabel('Pressure','FontSize',16)
% xlabel('Time','FontSize',16)
% subplot(313)
% plot(all_sd_tdc.')
% ax = gca;
% ax.FontSize = 16; 
% ylabel('Temperature','FontSize',16)
% xlabel('Time','FontSize',16)


% Plotting mean metrics of different objects
% close all
% 1 to display as subplots
% 0 to display as separate figures
subplots = 0;

% 1 to zoom between xlim([0 50])
zoom = 0;

% F0 High Frequency Fluid Vibrations
obj_name = {'Acrylic','Black foam','Car sponge','Flour sack','Kitchen sponge','Steel vase'};
plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#000000"];

if(subplots==0)
    figure();
else
    fig = figure('Position',  [100, 100, 1000, 900]);
    subplot(7,1,[1 2])
end
hold on
for j = 1:n_objects
    my_colour = plot_colours(j);
    plot(mean_objects{j}.vibrations,'Color',my_colour,'LineWidth',1.0)
end

if(zoom==1)
    xlim([0 50])
end
%ylim([1900 2150])
ax = gca;
ax.FontSize = 16;
title('Mean Vibration Time Series for each object','FontSize',16);
ylabel('High Frequency Fluid Vibrations','FontSize',16)
xlabel('Time','FontSize',16)
if(subplots==0)
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
% % saveas(gcf,'results/data_processing/vibrations_mean_timeseries.eps','epsc')
end

% F0 Low Frequency Fluid Pressure
if(subplots==0)
    figure();
else
    subplot(7,1,[3 4])
end
hold on
for j = 1:n_objects
    my_colour = plot_colours(j);
    plot(mean_objects{j}.pressure,'Color',my_colour,'LineWidth',1.0)
end

if(zoom==1)
    xlim([0 50])
end
ax = gca;
ax.FontSize = 16;
title('Mean Pressure Time Series for each object','FontSize',16);
ylabel('Low Frequency Fluid Pressure','FontSize',16)
xlabel('Time','FontSize',16)
if(subplots==0)
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
% % saveas(gcf,'results/data_processing/pressure_mean_timeseries.eps','epsc')
end

% F0 Core Temperature
if(subplots==0)
    figure();
else
    subplot(7,1,[5 6])
end
hold on
for j = 1:n_objects
    my_colour = plot_colours(j);
    plot(mean_objects{j}.temperature,'Color',my_colour,'LineWidth',1.0)
end

if(zoom==1)
    xlim([0 50])
end
ax = gca;
ax.FontSize = 16; 
title('Mean Temperature Time Series for each object','FontSize',16);
ylabel('Core Temperature','FontSize',16)
xlabel('Time','FontSize',16)
if(subplots==0)
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
% % saveas(gcf,'results/data_processing/temp_mean_timeseries.eps','epsc')
end

% if(subplots==1)
%     hSub = subplot(7,1,7);
%     hold on
%     for i = 1:n_objects
%         plot([1 nan],[1 nan],'Color',plot_colours(i),'LineWidth',1.0)
%     end
%     set(hSub, 'Visible', 'off');
%     legend(obj_name,'FontSize',16,'Location','southeastoutside')
% end

% Electrode impedance
figure();
count = 0;
for j = 1:n_objects
    count = count+1;
    %subplot(3,2,j)
    hold on
    my_colour = plot_colours(j);
    range = 1000*(j-1);
    for k = 1:19
        if k == 1
            plot([range:1:range+length(mean_objects{j}.electrodes(k,:))-1],mean_objects{j}.electrodes(k,:)...
            ,'Color',my_colour,'LineWidth',1.0);
        else
            h = plot([range:1:range+length(mean_objects{j}.electrodes(k,:))-1],mean_objects{j}.electrodes(k,:)...
            ,'Color',my_colour,'LineWidth',1.0);
            set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    ylim([2500 3800])
end
space = repelem('a',18);
for i = 1:(n_objects-1)
    xline(i*1000,'-k','LineWidth',1.0);
end
ax = gca;
ax.FontSize = 16; 
title('Mean Electrode Impedance Time Series for each object','FontSize',16);
ylabel('Electrode Impedance','FontSize',16)
xlabel('Time','FontSize',16)
legend(obj_name,'location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
XTicks = xticklabels;
newlab = cell(size(XTicks));
for i = 1:length(XTicks)
    newlab{i} = num2str(str2num(XTicks{i}) - floor(str2num(XTicks{i})/1000)*1000);
end
xticklabels(newlab);
% saveas(gcf,'results/data_processing/elec_mean_timeseries.eps','epsc')


% Plotting metrics of different objects and trials
close all
% F0 High Frequency Fluid Vibrations

% trials_to_plot = [1 2];
trials_to_plot = 1:1:10;
plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#000000"];

figure();
hold on
for i = 1:length(trials_to_plot)
    trial = trials_to_plot(i);
    for j = 1:n_objects
        my_colour = plot_colours(j);
        plot(objects{trial,j}.F0pac(2,:),'Color',my_colour,'LineWidth',1.0)
    end
end

% xlim([0 100])
% ylim([1900 2150])
ax = gca;
ax.FontSize = 16;
title('Vibrations Time Series for all trials and objects','FontSize',16);
ylabel('High Frequency Fluid Vibrations','FontSize',16)
xlabel('Time','FontSize',16)
xlim([0 1000])
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
% saveas(gcf,'results/data_processing/vibrations_trials_timeseries.eps','epsc')

% F0 Low Frequency Fluid Pressure
figure();
hold on
for i = 1:length(trials_to_plot)
    trial = trials_to_plot(i);
    for j = 1:n_objects
        my_colour = plot_colours(j);
        plot(objects{trial,j}.F0pdc,'Color',my_colour,'LineWidth',1.0)
    end
end

ax = gca;
ax.FontSize = 16;
title('Pressure Time Series for all trials and objects','FontSize',16);
ylabel('Low Frequency Fluid Pressure','FontSize',16)
xlabel('Time','FontSize',16)
xlim([0 1000])
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
% saveas(gcf,'results/data_processing/pressure_trials_timeseries.eps','epsc')

% F0 Core Temperature

figure();
hold on
for i = 1:length(trials_to_plot)
    trial = trials_to_plot(i);
    for j = 1:n_objects
        my_colour = plot_colours(j);
        plot(objects{trial,j}.F0tdc,'Color',my_colour,'LineWidth',1.0)
    end
end

ax = gca;
ax.FontSize = 16; 
title('Temperature Time Series for all trials and objects','FontSize',16);
ylabel('Core Temperature','FontSize',16)
xlabel('Time','FontSize',16)
xlim([0 1000])
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
% saveas(gcf,'results/data_processing/temp_trials_timeseries.eps','epsc')


% Electrode impedance

trials_to_plot = [1 2];
%trials_to_plot = 1:1:10;
plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#000000"];

figure();
count = 0;
for i = 1:length(trials_to_plot)
    trial = trials_to_plot(i);
    for j = 1:n_objects
        count = count+1;
        %subplot(3,2,j)
        hold on
        my_colour = plot_colours(j);
        range = 1000*(j-1);
        for k = 1:19
            if k == 1
                plot([range:1:range+length(mean_objects{j}.electrodes(k,:))-1],mean_objects{j}.electrodes(k,:)...
                ,'Color',my_colour,'LineWidth',1.0);
            else
                h = plot([range:1:range+length(mean_objects{j}.electrodes(k,:))-1],mean_objects{j}.electrodes(k,:)...
                ,'Color',my_colour,'LineWidth',1.0);
                set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
        ylim([2500 3800])
    end
end

for i = 1:(n_objects-1)
    xline(i*1000,'-k','LineWidth',1.0);
end

ax = gca;
ax.FontSize = 16; 
title('Electrode Impedance Time Series for all trials and objects','FontSize',16);
ylabel('Electrode Impedance','FontSize',16)
xlabel('Time','FontSize',16)
legend(obj_name,'FontSize',16,'Location','southeastoutside')
set(gcf, 'Position',  [100, 100, 900, 400])
XTicks = xticklabels;
newlab = cell(size(XTicks));
for i = 1:length(XTicks)
    newlab{i} = num2str(str2num(XTicks{i}) - floor(str2num(XTicks{i})/1000)*1000);
end
xticklabels(newlab);
% saveas(gcf,'results/data_processing/elec_trials_timeseries.eps','epsc')

% Section A - 2
% F0_PVT.mat and F0_Electrodes.mat

vibs = zeros(n_objects,trials);
pressures = zeros(n_objects,trials);
temps = zeros(n_objects,trials);
elecs = cell(n_objects,1);
for i = 1:n_objects
    vib = zeros(1,trials);
    pressure = zeros(1,trials);
    temp = zeros(1,trials);
    elec = zeros(19,trials);
    for j = 1:trials
        vib(j) = objects{j,i}.F0pac(2,time_step);
        pressure(j) = objects{j,i}.F0pdc(time_step);
        temp(j) = objects{j,i}.F0tdc(time_step);
        elec(:,j) = objects{j,i}.F0Electrodes(:,time_step);
    end
    vibs(i,:) = vib;
    pressures(i,:) = pressure;
    temps(i,:) = temp;
    elecs{i} = elec;
end

acrylic_pvt = struct('vibrations',vibs(1,:),'pressure',pressures(1,:),...
    'temp',temps(1,:));
black_foam_pvt = struct('vibrations',vibs(2,:),'pressure',pressures(2,:),...
    'temp',temps(2,:));
car_sponge_pvt = struct('vibrations',vibs(3,:),'pressure',pressures(3,:),...
    'temp',temps(3,:));
flour_sack_pvt = struct('vibrations',vibs(4,:),'pressure',pressures(4,:),...
    'temp',temps(4,:));
kitchen_sponge_pvt = struct('vibrations',vibs(5,:),'pressure',pressures(5,:),...
    'temp',temps(5,:));
steel_vase_pvt = struct('vibrations',vibs(6,:),'pressure',pressures(6,:),...
    'temp',temps(6,:));

acrylic_electrodes = elecs{1};
black_foam_electrodes = elecs{2};
car_sponge_electrodes = elecs{3};
flour_sack_electrodes = elecs{4};
kitchen_sponge_electrodes = elecs{5};
steel_vase_electrodes = elecs{6};

% save('data-matrices/F0_PVT.mat','acrylic_pvt','black_foam_pvt','car_sponge_pvt',...
%     'flour_sack_pvt','kitchen_sponge_pvt','steel_vase_pvt');
% 
% save('data-matrices/F0_Electrodes.mat','acrylic_electrodes','black_foam_electrodes',...
%     'car_sponge_electrodes','flour_sack_electrodes',...
%     'kitchen_sponge_electrodes','steel_vase_electrodes');

% Section A - 3
% Scatter plot
obj_name = {'Acrylic','Black foam','Car sponge','Flour sack','Kitchen sponge','Steel vase'};
trials = 10;
n_objects = 6;

load('data-matrices/F0_PVT.mat')
load('data-matrices/F0_Electrodes.mat')

all_pvt = {acrylic_pvt,black_foam_pvt,car_sponge_pvt,...
    flour_sack_pvt,kitchen_sponge_pvt,steel_vase_pvt};

all_electrodes = {acrylic_electrodes,black_foam_electrodes,...
    car_sponge_electrodes,flour_sack_electrodes,...
    kitchen_sponge_electrodes,steel_vase_electrodes};

plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#000000"];

figure();
subplot(1,5,[1 4])
for i = 1:length(all_pvt)
    scatter3(all_pvt{i}.pressure,all_pvt{i}.vibrations,all_pvt{i}.temp,75,...
        'o','filled','MarkerFaceColor',plot_colours(i))
    hold on;
end
ax = gca;
ax.FontSize = 16;
title('3D Scatter Plot of the PVT data for all trials and objects','FontSize',16);
xlabel('Pressure','FontSize',16)
ylabel('Vibrations','FontSize',16)
zlabel('Temperature','FontSize',16)

hSub = subplot(1,5,5);
hold on
for i = 1:n_objects
    plot([1 nan],[1 nan],'Color',plot_colours(i),'LineWidth',1.0)
end
set(hSub, 'Visible', 'off');
legend(obj_name,'FontSize',16,'Location','eastoutside')
set(gcf, 'Position',  [100, 0, 1000, 600])
% saveas(gcf,'results/data_processing/pvt_3d_plot.eps','epsc')
