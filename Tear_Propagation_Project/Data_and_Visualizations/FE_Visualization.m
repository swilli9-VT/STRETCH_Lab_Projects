%% Reset
clear;
close all;
clc;
%% Retrieve mesh info

job_directory = "D:\Users\Will\Tear_Propagation_Project\";
job_name = "HGO_Tear_Propagation";

%Define initial nodal positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(job_directory+job_name+".inp");
item1 = '*Node';
nextItem(fileID,item1);
coord = textscan(fileID, '%u %f %f %f %*s','delimiter',',');
fclose(fileID);
    
X = coord{1,2};
Y = coord{1,3};
Z = coord{1,4};

coordinates = [X Y Z];
num_dof = length(coordinates)*3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Get nodal connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(job_directory+job_name+".inp");
item1 = '*Element, type=';
nextItem(fileID,item1);
connect = textscan(fileID, '%u %u %u %u %u %u %u %u %u','delimiter',',\n');
fclose(fileID);

nodes = zeros(length(connect{1}),length(connect)-1);

for i = 1:size(nodes,2)
    nodes(:,i) = connect{1,i+1};
end

%% Load data from CSVs

Data_Snapshots = table2array(readtable('HGO_Tear_Propagation_Snapshots.csv'));
Data_Status = table2array(readtable('HGO_Tear_Propagation_Status.csv'));
% Data_Beta = table2array(readtable('HGO_Tear_Propagation_Beta.csv')); 
% Data_Alpha = table2array(readtable('HGO_Tear_Propagation_E_Alpha.csv'));
% Data_Damage = table2array(readtable('HGO_Tear_Propagation_Damage.csv'));
% Data_Strain = table2array(readtable('HGO_Tear_Propagation_Strain.csv'));

%% Sort data by mu parameter case

beta_d = [16.75, 40.78];
beta_m = [26.11, 18.85];
beta_p = [69.21, 63.2];

mu_cases = [];

mu_times = cell(1,combination);
mu_snapshots = cell(1,combination);
mu_status = cell(1,combination);
% mu_beta = cell(1,combination);
% mu_dam = cell(1,combination);
% mu_alpha = cell(1,combination);
% mu_strain = cell(1,combination);

c=1;
for b1 = 1:length(beta_d)
    for b2 = 1:length(beta_m)
        for b3 = 1:length(beta_p)
            snap_mu = Data_Snapshots((Data_Snapshots(:,2)==beta_d(b1) &...
                 Data_Snapshots(:,3)==beta_m(b2) &...
                 Data_Snapshots(:,4)==beta_p(b3)), :);

            stat_mu = Data_Status((Data_Status(:,1)==beta_d(b1) &...
                 Data_Status(:,2)==beta_m(b2) &...
                 Data_Status(:,3)==beta_p(b3)), :);


            mu_times{1,c} = snap_mu(:,5)';
            mu_snapshots{1,c} = snap_mu(:,6:end)';
            mu_status{1,c} = stat_mu(:,5:end)';

            mu_cases = [mu_cases; beta_d(b1), beta_m(b2), beta_p(b3)];
            c = c+1;
        end
    end
end

clearvars -except coordinates nodes mu_times mu_snapshots mu_status mu_cases

%% Make still-image plots of data
close all;

combination = length(mu_cases);
num_pressures_to_plot = 1;
factor = 1; %set scale factor

fig=1;
for i = 1:combination
    pressures_index = zeros(1,num_pressures_to_plot);
    intervals = [1];%linspace(0,1,num_pressures_to_plot);
    for t = 1:length(intervals)
        [~, pressures_index(t)] = min(abs(mu_times{i}-intervals(t)));
    end
    
    for p = 1:num_pressures_to_plot
        U_comp = zeros(length(mu_snapshots{i})/3, 3);
        x=1;
        for q = 1:3:length(mu_snapshots{i})
            U_comp(x,1) = mu_snapshots{i}(q,pressures_index(p));
            U_comp(x,2) = mu_snapshots{i}(q+1,pressures_index(p));
            U_comp(x,3) = mu_snapshots{i}(q+2,pressures_index(p));
            x = x+1;
        end
        components = [U_comp(:,1) U_comp(:,2) U_comp(:,3)];
        magnitude = sqrt(U_comp(:,1).^2+U_comp(:,2).^2+U_comp(:,3).^2);
        el_active = logical(mu_status{i}(:,pressures_index(p)));

        figure(fig)
        %subplot(1,size(deformSets,2),i)
        PlotFieldonDefoMesh(coordinates,nodes(el_active,:),...
                        factor,components,magnitude,...
                        [0 4.1], "linear", '$||u||$ (mm)');
        axis equal
        hold on
        quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1])
        %title('$t='+string(mu_times{i}(pressures_index(p)))+'$',interpreter="latex")
        title("$\mu_{"+string(i)+"} = \{"+string(mu_cases(i,1))+...
                ", "+string(mu_cases(i,2))+", "+string(mu_cases(i,3))+"\}$", interpreter="latex")
        xlabel('x (mm)')
        ylabel('y (mm)')
        zlabel('z (mm)')
        axis([-8 8 -8 8 -1 15])
        xticks([-8 -6 -4 -2 0 2 4 6 8])
        yticks([-8 -6 -4 -2 0 2 4 6 8])
        zticks([0 2 4 6 8 10 12 14])
        set(gca, 'xdir','reverse','ydir','reverse')
        set(gca, 'fontsize',14,'FontName', 'Times','XColor', [0 0 0],'YColor', [0 0 0])
        axes.color = 'black';
        %colorbar(gca, 'off')
        %view([0 0 1])
        hold off

        fig = fig+1;

    end    
end

%% Make still-image plot of undeformed mesh
close all;

clc;
a_o = [4.71 1.95];
a_i = [4.29 2.71];

b_o = flip(a_o);
b_i = flip(a_i);

A_ixB_i = (a_i(1)*b_i(2) - a_i(2)*b_i(1)); % always positive
A_oxB_o = (a_o(1)*b_o(2) - a_o(2)*b_o(1)); % always positive

ind_away_from_tear = [];
ind_near_tear = [];
for i=1:length(coordinates)
    c = coordinates(i,1:2);
    
    A_ixC = (a_i(1)*c(2) - a_i(2)*c(1));
    CxB_i = (c(1)*b_i(2) - c(2)*b_i(1));
    CxA_i = (c(1)*a_i(2) - c(2)*a_i(1));

    A_oxC = (a_o(1)*c(2) - a_o(2)*c(1));
    CxB_o = (c(1)*b_o(2) - c(2)*b_o(1));
    CxA_o = (c(1)*a_o(2) - c(2)*a_o(1));

    if ~(CxB_i * CxA_i <= 0 && A_ixB_i * A_ixC >=0)
        ind_away_from_tear = [ind_away_from_tear i];
    end

    if (CxB_o * CxA_o <= 0 && A_oxB_o * A_oxC >=0)
        ind_near_tear = [ind_near_tear i];
    end
end

nodes_away = [];
nodes_near = [];
magnitude = zeros(length(coordinates),1);
for i=1:length(nodes)
    [found_a, ~] = ismember(ind_away_from_tear, nodes(i,:));
    [found_n, ~] = ismember(ind_near_tear, nodes(i,:));

    flag_a = sum(found_a) == 8;
    flag_n = sum(found_n) == 8;

    if (flag_a && ~flag_n)
        nodes_away = [nodes_away; nodes(i,:)];
        %magnitude(nodes(i,:),:) = 1;
    
    elseif (flag_n && ~flag_a)
        nodes_near = [nodes_near; nodes(i,:)];
        %magnitude(nodes(i,:),:) = 2;

    elseif (flag_a && flag_n)
        nodes_away = [nodes_away; nodes(i,:)];
        nodes_near = [nodes_near; nodes(i,:)];
    end
end

overlap = intersect(nodes_away, nodes_near, 'rows');

nodes_overlap = nodes(ismember(nodes,overlap,'rows'),:);

% 
% overlap = []
% for i=1:length(nodes_away)
%     for j=1:length(nodes_near)
%         if nodes_away(i,:)==nodes_near
%     end
% end

factor = 1; %set scale factor

fig=1;
i=1;

U_comp = zeros(length(mu_snapshots{i})/3, 3);
components = [U_comp(:,1) U_comp(:,2) U_comp(:,3)];
el_active = logical(mu_status{i}(:,1));

figure(fig)
%subplot(1,size(deformSets,2),i)
PlotFieldonDefoMesh(coordinates,nodes_overlap,...
                factor,components,magnitude,...
                [0 2.0], "linear", '$||u||$ (mm)');
axis equal
hold on
quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1])
%title(sublabels(i))
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
axis([-8 8 -8 8 -1 15])
xticks([-8 -6 -4 -2 0 2 4 6 8])
yticks([-8 -6 -4 -2 0 2 4 6 8])
zticks([0 2 4 6 8 10 12 14])
set(gca, 'xdir','reverse','ydir','reverse')
set(gca, 'fontsize',14,'FontName', 'Times','XColor', [0 0 0],'YColor', [0 0 0])
axes.color = 'black';
%colorbar(gca, 'off')
%view([0 0 1])
hold off


%% Animate data
close all;

factor = 1; %set scale factor

clearvars vidObj

for i = 1:combination
    num_steps = size(mu_snapshots{i},2);

    avi_str = 'HGO_Tear_Propagation'+string(i)+'.avi';

%     M(num_steps) = struct('cdata',[],'colormap',[]);
    vidObj = VideoWriter(avi_str);
    vidObj.open()

    h = figure(i);
    h.Visible = 'off';

    hold on
    xlabel('x (mm)')
    ylabel('y (mm)')
    zlabel('z (mm)')
    axis equal
    axis([-8 8 -8 8 -1 15])
    xticks([-8 -6 -4 -2 0 2 4 6 8])
    yticks([-8 -6 -4 -2 0 2 4 6 8])
    zticks([0 2 4 6 8 10 12 14])
    set(gca, 'xdir','reverse','ydir','reverse')
    set(gca, 'fontsize',14,'FontName', 'Times','XColor', [0 0 0],'YColor', [0 0 0])
    axes.color = 'black';
    %colorbar(gca, 'off')
    
    for p = 1:10
        U_comp = zeros(length(mu_snapshots{i})/3, 3);
        x=1;
        for q = 1:3:length(mu_snapshots{i})
            U_comp(x,1) = mu_snapshots{i}(q,p);
            U_comp(x,2) = mu_snapshots{i}(q+1,p);
            U_comp(x,3) = mu_snapshots{i}(q+2,p);
            x = x+1;
        end
        components = [U_comp(:,1) U_comp(:,2) U_comp(:,3)];
        magnitude = sqrt(U_comp(:,1).^2+U_comp(:,2).^2+U_comp(:,3).^2);
        el_active = logical(mu_status{i}(:,p));
        
        title('$t='+string(mu_times{i}(p))+'$',interpreter="latex")
        quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1])
        PlotFieldonDefoMesh(coordinates,nodes(el_active,:),...
                        factor,components,magnitude,...
                        [0 4.1], "linear", '$||u||$ (mm)');

        if p == 1
            saveas(h,'PlayFrame_'+string(i)+'.png')
        end
        currFrame = getframe(h);
        writeVideo(vidObj,currFrame);
        cla;

        percent_complete = 100*(p/num_steps);
        
        clc;
        disp("Movie "+string(i)+" progress: "+string(percent_complete)+"%")
    end 

    hold off
    for f=1:60 % write another 2 seconds of the last frame onto the end of the video
        writeVideo(vidObj,currFrame);
    end
    close(vidObj);
end


%% Test Section

test = [1 2 3; 3 2 1; 2 1 3; 3 1 2];
check = [1 2 3];

ismember(test, check, 'rows')

%disp(where)

% fullsetcount = accumarray(where(found)', 1, [numel(basesetvalues), 1]);
% diffcount = basesetcount - fullsetcount;
% notenough = diffcount > 0;
% out = table(basesetvalues(notenough).', diffcount(notenough), 'VariableNames', {'Value', 'MissingCount'})