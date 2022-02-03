% Section C

clear all; close all; clc;

trials = 10;
n_objects = 6;

load('data-matrices/F0_PVT.mat')
load('data-matrices/F0_Electrodes.mat')

all_objects = cell(1,6);
all_objects{1} = acrylic_pvt;
all_objects{2} = black_foam_pvt;
all_objects{3} = car_sponge_pvt;
all_objects{4} = flour_sack_pvt;
all_objects{5} = kitchen_sponge_pvt;
all_objects{6} = steel_vase_pvt;

object1 = 2;
object2 = 3;
obj_name = {'Black foam','Car sponge'};
obj_label = 'bc';

% object1 = 4;
% object2 = 5;
% obj_name = {'Flour sack','Kitchen sponge'};
% obj_label = 'fk';

plot_colours = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#000000"];
plot_colours2 = ["0072BD","#ff0000","#ed8020","#7b00ff","#8beb0c","#000000"];
plot_colours = [plot_colours(object1) plot_colours(object2)];
plot_colours2 = [plot_colours2(object1) plot_colours2(object2)];

my_objects = cell(3,2);
my_objects{1,1} = all_objects{object1}.pressure;
my_objects{2,1} = all_objects{object1}.vibrations;
my_objects{3,1} = all_objects{object1}.temp;
my_objects{1,2} = all_objects{object2}.pressure;
my_objects{2,2} = all_objects{object2}.vibrations;
my_objects{3,2} = all_objects{object2}.temp;

N = length(acrylic_pvt.pressure);

% 2 features comparison
compare = [1 2; 1 3; 3 2];  % which pair of features
y = [zeros(1,N) ones(1,N)];

data = cell(1,3);
models = cell(1,3);

% gradient_offset = [0.5 -0.55 0.15]; % pv = 0.5, pt = 0.15, tv = -0.55
for i = 1:length(compare)
    temp = zeros(2*N,2);
    temp(1:N,1) = my_objects{compare(i,1),1};
    temp(1:N,2) = my_objects{compare(i,2),1};
    
    temp(N+1:end,1)= my_objects{compare(i,1),2};
    temp(N+1:end,2)= my_objects{compare(i,2),2};
    
    temp = zscore(temp);
    
    mean1 = mean(temp(1:N,:),1);
    mean2 = mean(temp(N+1:end,:),1);
    mean_overall = mean(temp,1);
    
    temp_mat = temp(1:N,:);
    Sw = temp_mat - mean1;
    Sw = Sw' * Sw;
    
    temp_mat = temp(N+1:end,:);
    Sw2 = temp_mat - mean2;
    Sw2 = Sw2' * Sw2;
    Sw = Sw + Sw2;
    
    temp_mat = mean1 - mean2;
    Sb = temp_mat' * temp_mat;
    
    temp_mat = inv(Sw)*Sb;
    [V,D] = eig(temp_mat);
    D = sum(D,1);
    
    max_idx = find(D==max(D));
    my_vector = V(:,max_idx);
    scaling = 100;
    my_vector = scaling .* my_vector;
    sep_line = [-my_vector(2) my_vector(1)];
    
    data{i} = temp;

    models{i} = fitcdiscr(temp,y);
    
    x1plot = linspace(1.01*min(temp(:,1)), 1.01*max(temp(:,1)), 500)';
    x2plot = linspace(1.01*min(temp(:,2)), 1.01*max(temp(:,2)), 500)';
    [X1, X2] = meshgrid(x1plot, x2plot);
    vals = zeros(size(X1));
    for j = 1:size(X1, 2)
       this_X = [X1(:,j), X2(:,j)];
       vals(:,j) = predict(models{i}, this_X);
    end
    
    figure();
    hold on
    plot(temp(y==0,2),temp(y==0,1),'+','LineWidth',1.5,'Color',plot_colours(1));
    plot(temp(y==1,2),temp(y==1,1),'*','LineWidth',1.5,'Color',plot_colours(2));

    scatter(mean1(2),mean1(1),75,...
        'd','filled','MarkerFaceColor',plot_colours2(1))
    
    scatter(mean2(2),mean2(1),75,...
        'd','filled','MarkerFaceColor',plot_colours2(2))
    
    scatter(mean_overall(2),mean_overall(1),75,...
        'o','filled','MarkerFaceColor','black')
    
    plot([-my_vector(2) my_vector(2)],[-my_vector(1) my_vector(1)],'g-','LineWidth',1.5)
 
    plot([-sep_line(2) sep_line(2)],[-sep_line(1) sep_line(1)],'b-.','LineWidth',1.5)
    
    legend(obj_name,'FontSize',16,'Location','southwest')
    ax = gca;
    ax.FontSize = 16; 
    
    text = cell(1,2);
    feature = cell(1,2);
    for j = 1:length(compare(i,:))
        if(compare(i,j)==1)
            text{j} = "Pressure";
            feature{j} = "p";
        elseif(compare(i,j)==2)
            text{j} = "Vibrations";
            feature{j} = "v";
        else
            text{j} = "Temperature";
            feature{j} = "t";
        end
    end
    
    title(sprintf('LD function (Green) between %s and %s',text{1},text{2}),'FontSize',14)
    ylabel(text(1),'FontSize',16)
    xlabel(text(2),'FontSize',16)
    ylim([1.1*min(temp(:,1)) 1.1*max(temp(:,1))])
    xlim([1.1*min(temp(:,2)) 1.1*max(temp(:,2))])
    % saveas(gcf,sprintf('results/lda/%s_lda_2features_%s%s_2d.eps',obj_label,feature{1},feature{2}),'epsc')
end


% 1D number line projection onto LDA (Green) for 2 features comparison 
figure();
for i = 1:length(compare)
    temp = zeros(2*N,2);
    temp(1:N,1) = my_objects{compare(i,1),1};
    temp(1:N,2) = my_objects{compare(i,2),1};

    temp(N+1:end,1)= my_objects{compare(i,1),2};
    temp(N+1:end,2)= my_objects{compare(i,2),2};
    
    temp = zscore(temp);
    
    mean1 = mean(temp(1:N,:),1);
    mean2 = mean(temp(N+1:end,:),1);
    mean_overall = mean(temp,1);
    
    temp_mat = temp(1:N,:);
    Sw = temp_mat - mean1;
    Sw = Sw' * Sw;
    
    temp_mat = temp(N+1:end,:);
    Sw2 = temp_mat - mean2;
    Sw2 = Sw2' * Sw2;
    Sw = Sw + Sw2;
    
    temp_mat = mean1 - mean2;
    Sb = temp_mat' * temp_mat;
    
    temp_mat = inv(Sw)*Sb;
    [V,D] = eig(temp_mat);
    D = sum(D,1);
    
    max_idx = find(D==max(D));
    my_vector = V(:,max_idx);
    scaling = 100;
    my_vector = scaling .* my_vector;
   
    % Projection onto LD green (1D number line)
    proj_points = zeros(length(temp),2);
    proj_mean1 = proj1D([-my_vector(1) -my_vector(2) ; my_vector(1) my_vector(2)],mean1);
    proj_mean2 = proj1D([-my_vector(1) -my_vector(2) ; my_vector(1) my_vector(2)],mean2);
    for k = 1:length(temp)
        proj_points(k,:) = proj1D([-my_vector(1) -my_vector(2) ; my_vector(1) my_vector(2)],temp(k,:));
    end
    
    hold on
    h1 = xline(0,'k--','linewidth',1.1);
    h2 = scatter(proj_mean1(1),4-i,100,'d','filled','MarkerFaceColor',plot_colours2(1));
    h3 = scatter(proj_mean2(1),4-i,100,'d','filled','MarkerFaceColor',plot_colours2(2));
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(proj_points(y==0),4*ones(N,1)-i,'+','LineWidth',1.5,'Color',plot_colours(1));
    plot(proj_points(y==1),4*ones(N,1)-i,'*','LineWidth',1.5,'Color',plot_colours(2));
end
title('Projection of 2-features PVT data onto the LD function','FontSize',14);
ylim([0 4])
yticks([1 2 3])
yticklabels({'TV','PT','PV'});
ylabel('Feature Space')
xlabel('Score','FontSize',16);
ax = gca;
ax.YAxis.FontSize = 16;
legend show
legend(obj_name,'FontSize',16,'Location','southwest')
% saveas(gcf,sprintf('results/lda/%s_lda_2features_1d.eps',obj_label),'epsc')

% 3 features comparison

data = zeros(2*N,3);
data(1:N,1) = my_objects{1,1};
data(1:N,2) = my_objects{2,1};
data(1:N,3) = my_objects{3,1};

data(N+1:end,1)= my_objects{1,2};
data(N+1:end,2)= my_objects{2,2};
data(N+1:end,3)= my_objects{3,2};

data = zscore(data);

mean1 = mean(data(1:N,:),1);
mean2 = mean(data(N+1:end,:),1);
mean_overall = mean(data,1);

temp_mat = data(1:N,:);
Sw = temp_mat - mean1;
Sw = Sw' * Sw;

temp_mat = data(N+1:end,:);
Sw2 = temp_mat - mean2;
Sw2 = Sw2' * Sw2;
Sw = Sw + Sw2;

temp_mat = mean1 - mean2;
Sb = temp_mat' * temp_mat;

temp_mat = inv(Sw)*Sb;
[V,D] = eig(temp_mat);
V = real(V);
D = real(D);
D = sum(D,1);

max_values = maxk(D,2);
max_idx = find(D==max_values(1),1);
max2_idx = find(D==max_values(2),1);
my_vector = V(:,max_idx);
my_vector2 = V(:,max2_idx);

scaling = 100;
my_vector = scaling .* my_vector;
my_vector2 = scaling .* my_vector2;
sep_line = [-my_vector(2) my_vector(1)];

y = [zeros(1,N) ones(1,N)];
% model3 = fitcdiscr(data,y);

% x1plot = linspace(1.1*min(data(:,1)), 1.1*max(data(:,1)), 10)';
% x2plot = linspace(1.1*min(data(:,2)), 1.1*max(data(:,2)), 10)';
% x3plot = linspace(1.1*min(data(:,3)), 1.1*max(data(:,3)), 10)';
% 
% [X1,X2,X3] = meshgrid(x1plot,x2plot,x3plot);
% vals = zeros(size(X1));
% 
% for i = 1:size(X1,1)
%     for j = 1:size(X2,1)
%         for k = 1:size(X3,1)
%             this_X = [X1(i,j,k) X2(i,j,k) X3(i,j,k)];
%             vals(i,j,k) = predict(model3, this_X);
%         end
%     end
% end

figure();
subplot(1,5,[1 4])
hold on;

% x_axis = [-my_vector(1) -my_vector2(1) my_vector(1) my_vector2(1)];
% y_axis = [-my_vector(2) -my_vector2(2) my_vector(2) my_vector2(2)];
% z_axis = [-my_vector(3) -my_vector2(3) my_vector(3) my_vector2(3)];
% p = patch(isosurface(X1,X2,X3,vals));
% p = patch(x_axis,y_axis,z_axis);

v = [-my_vector2(1) -my_vector2(2) -my_vector2(3);...
    my_vector(1) my_vector(2) my_vector(3);...
    my_vector2(1) my_vector2(2) my_vector2(3);...
    -my_vector(1) -my_vector(2) -my_vector(3)];

f = [1 2 3 4];

p = patch('Faces',f,'Vertices',v);
p.FaceColor = '#CDCDCD';
p.EdgeColor = 'none';
p.FaceAlpha = 0.7;
p.EdgeAlpha = 1;

scatter3(data(y==0,1),data(y==0,2),data(y==0,3),75,...
    '+','MarkerEdgeColor',plot_colours(1),'LineWidth',1.1);
scatter3(data(y==1,1),data(y==1,2),data(y==1,3),75,...
    '*','MarkerEdgeColor',plot_colours(2),'LineWidth',1.1);

h1 = scatter3(mean1(2),mean1(1),mean1(3),75,...
       'd','filled','MarkerFaceColor',plot_colours2(1));
h2 = scatter3(mean2(2),mean2(1),mean2(3),75,...
        'd','filled','MarkerFaceColor',plot_colours2(2));
h3 = scatter3(mean_overall(2),mean_overall(1),mean_overall(3),75,...
        'o','filled','MarkerFaceColor','black');
    
h4 = plot3([-my_vector(1) my_vector(1)],[-my_vector(2) my_vector(2)],...
    [-my_vector(3) my_vector(3)],'g-','LineWidth',1.5);
h5 = plot3([-my_vector2(1) my_vector2(1)],[-my_vector2(2) my_vector2(2)],...
    [-my_vector2(3) my_vector2(3)],'b-','LineWidth',1.5);

set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

grid on;
view(3);
%view(90,0); % Z-Y plane
%view(2);    % Y-X plane
%view(180,0); % Z-X plane
%view(-90.214657585129416,-0.133131549117492);
ax = gca;
ax.FontSize = 16; 
title('LD functions for the 3D PVT data','FontSize',18)
xlabel('Pressure','FontSize',16)
ylabel('Vibrations','FontSize',16)
zlabel('Temperature','FontSize',16)
xlim([1.1*min(data(:,1)) 1.1*max(data(:,1))])
ylim([1.1*min(data(:,2)) 1.1*max(data(:,2))])
zlim([1.1*min(data(:,3)) 1.1*max(data(:,3))])
set(gcf, 'Position',  [100, 0, 1300, 900])
hSub = subplot(1,5,5);
hold on
plot([1 nan],[1 nan],'Color',plot_colours(1),'LineWidth',3)
plot([1 nan],[1 nan],'Color',plot_colours(2),'LineWidth',3)
set(hSub, 'Visible', 'off');
legend(obj_name,'FontSize',16,'Location','east')
% saveas(gcf,sprintf('results/lda/%s_lda_3features_3d.eps',obj_label),'epsc')

proj_points = zeros(size(data));  % project onto green axis
proj_points2 = zeros(size(data)); % project onto blue axis
for k = 1:length(data)
    proj_points(k,:) = proj2D([-my_vector(1) -my_vector(2) -my_vector(3);...
        my_vector(1) my_vector(2) my_vector(3)],data(k,:));
    proj_points2(k,:) = proj2D([-my_vector2(1) -my_vector2(2) -my_vector2(3);...
        my_vector2(1) my_vector2(2) my_vector2(3)],data(k,:));
end
proj_mean1_g = proj2D([-my_vector(1) -my_vector(2) -my_vector(3);...
        my_vector(1) my_vector(2) my_vector(3)],mean1);
proj_mean1_b  = proj2D([-my_vector2(1) -my_vector2(2) -my_vector2(3);...
        my_vector2(1) my_vector2(2) my_vector2(3)],mean1);
proj_mean2_g = proj2D([-my_vector(1) -my_vector(2) -my_vector(3);...
        my_vector(1) my_vector(2) my_vector(3)],mean2);
proj_mean2_b  = proj2D([-my_vector2(1) -my_vector2(2) -my_vector2(3);...
        my_vector2(1) my_vector2(2) my_vector2(3)],mean2);


% convert projected points to 1d coordinate
proj_points = proj_points*ones(3,1);
proj_points2 = proj_points2*ones(3,1);
proj_mean1_g = proj_mean1_g*ones(3,1);
proj_mean1_b = proj_mean1_b*ones(3,1);
proj_mean2_g = proj_mean2_g*ones(3,1);
proj_mean2_b = proj_mean2_b*ones(3,1);

% 1D projection
figure();
hold on
h6 = xline(0,'k--','LineWidth',1.1);
set(get(get(h6,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
scatter(proj_mean1_g,1,75,'d','filled','MarkerFaceColor',plot_colours2(1));
scatter(proj_mean2_g,1,75,'d','filled','MarkerFaceColor',plot_colours2(2));
plot(proj_points(y==0),ones(N,1),'+','MarkerSize',8,'LineWidth',1.5,'Color',plot_colours(1));
plot(proj_points(y==1),ones(N,1),'*','MarkerSize',8,'LineWidth',1.5,'Color',plot_colours(2));

scatter(proj_mean1_b,0,100,'d','filled','MarkerFaceColor',plot_colours2(1));
scatter(proj_mean2_b,0,100,'d','filled','MarkerFaceColor',plot_colours2(2));
plot(proj_points2(y==0),zeros(N,1),'+','MarkerSize',8,'LineWidth',1.5,'Color',plot_colours(1));
plot(proj_points2(y==1),zeros(N,1),'*','MarkerSize',8,'LineWidth',1.5,'Color',plot_colours(2));
ylim([-1 2])
yticks([0 1])
ylabel('LD Functions')
yticklabels({'Blue','Green'})
xlabel('Score','FontSize',16)
ax = gca;
ax.YAxis.FontSize = 16;
title('Projection of 3-features PVT data onto the LD functions','FontSize',14)
legend(obj_name,'FontSize',16,'Location','southeast')
% saveas(gcf,sprintf('results/lda/%s_lda_3features_1d.eps',obj_label),'epsc')

% 2D projection
figure()
hold on
scatter(proj_mean1_g,proj_mean1_b,100,'d','filled','MarkerFaceColor',plot_colours2(1));
scatter(proj_mean2_g,proj_mean2_b,100,'d','filled','MarkerFaceColor',plot_colours2(2));
plot(proj_points(y==0),proj_points2(y==0),'+','MarkerSize',8,'LineWidth',1.5,'Color',plot_colours(1));
plot(proj_points(y==1),proj_points2(y==1),'*','MarkerSize',8,'LineWidth',1.5,'Color',plot_colours(2));
xlabel('LD (Green)','FontSize',16);
ylabel('LD (Blue)','FontSize',16);
title('Projection of 3-features PVT data onto the plane defined by LD functions','FontSize',14)
legend(obj_name,'FontSize',16,'Location','southeast')
% saveas(gcf,sprintf('results/lda/%s_lda_3features_2d.eps',obj_label),'epsc')

%
function [ProjPoint] = proj1D(vector, q)
p0 = vector(1,:);
p1 = vector(2,:);
a = [-q(1)*(p1(1)-p0(1)) - q(2)*(p1(2)-p0(2)); ...
    -p0(2)*(p1(1)-p0(1)) + p0(1)*(p1(2)-p0(2))]; 
b = [p1(1) - p0(1), p1(2) - p0(2);...
    p0(2) - p1(2), p1(1) - p0(1)];
ProjPoint = -(b\a);
end

function [ProjPoint] = proj2D(vector, p)
a = vector(1,:);
b = vector(2,:);
ap = p-a;
ab = b-a;
ProjPoint = a + dot(ap,ab)/dot(ab,ab) * ab;
end


