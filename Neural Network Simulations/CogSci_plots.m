% author: Sebastian Musslick

%% PLOT MIS VALIDATION
%% extract relevant data

MSE_data = nan(length(batch_log(1,1).multiPerformance_mean), 2, NPathways);

for rep = 1:length(batch_log(1,1).multiPerformance_mean)
    
    for capIdx = 1:size(batch_log(1,1).multiPerformance_mean{rep},2)
       
        MSE_data(rep, 1, batch_log(1,1).multiPerformance_mean{rep}(1, capIdx)) = batch_log(1,1).multiPerformance_mean{rep}(2, capIdx); % independent sets
        MSE_data(rep, 2, batch_log(1,1).multiPerformance_mean{rep}(1, capIdx)) = batch_log(1,1).multiPerformance_mean{rep}(3, capIdx); % dependent sets
        
    end
    
end

% compress data

MSE_data_sem = zeros(3, NPathways);
MSE_data_mean = zeros(3, NPathways);

MSE_data_sem(1,:) = 1:NPathways;
MSE_data_mean(1,:) = 1:NPathways;

MSE_data_mean(2:3,:) = squeeze(nanmean(MSE_data, 1));

for cap = 1:size(MSE_data_mean, 2)
    
    MSE_data_sem(2, cap) = nanstd(MSE_data(:,1,cap))/sqrt(sum(~isnan(MSE_data(:,1,cap))));
    MSE_data_sem(3, cap) = nanstd(MSE_data(:,2,cap))/sqrt(sum(~isnan(MSE_data(:,2,cap))));
    
end

% clean for NaN

nonNanIdx = find(sum(isnan(MSE_data_mean)) == 0);

MSE_data_mean = MSE_data_mean(:,nonNanIdx);
MSE_data_sem = MSE_data_sem(:,nonNanIdx);

%% plot data

colorA = [0 128 0]/255;
colorB = [0 0 128]/255;
colorC = [102 0 102]/255;
colorD = [128 0 0]/255;

labelSize = 18;
lineWidth = 5;
fontName = 'Helvetica';

color = [50 50 50;150 150 150]/255; % nice color: 0 51 102
color = [0 0 100; 100 0 0]/255;

bardata_mean = MSE_data_mean(2:3,1:3)';
            
bardata_sem = MSE_data_sem(2:3,1:3)';

%bar_handle = bar(1:4,bardata_mean, 'FaceColor',color); hold on;
xlim([0.5 size(bardata_mean,1)+0.5]);
bar_handle = errorbar_groups(bardata_mean', bardata_sem','bar_colors',color,'errorbar_width',0.5,'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',1.5});
set(gca,'XTickLabel',{' ', ' ',' ',' '},'FontSize', labelSize, 'FontName', fontName);
ylabel('MSE','FontSize', labelSize, 'FontName', fontName);
xlabel('Total Number Of Tasks Performed','FontSize', labelSize);
%title({'multitasking performance'},'FontSize', labelSize, 'FontName', fontName);
l = legend('Independent','Non-Independent','Location','northwest');
set(l, 'FontName', fontName);
%ylim([0 0.07]);
set(gcf, 'Position', [800 800 380 250])
set(gca, 'XTickLabel', MSE_data_mean(1,:));
set(gcf, 'Color', [1 1 1])

%% corresponding t-test


for cap = 2:size(MSE_data,3)
   
    MSE_independent = MSE_data(~isnan(MSE_data(:,1,cap)),1,cap);
    MSE_dependent = MSE_data(~isnan(MSE_data(:,2,cap)),2,cap);
    
    [H,P,CI,STATS] = ttest(MSE_dependent-MSE_independent, 0, 0.05, 'right');
    disp(['For cap = ' num2str(cap) ': t = ' num2str(STATS.tstat) ', d = ' num2str(STATS.df) ', sd = ' num2str(STATS.sd), ', mean = ' num2str(mean(MSE_dependent-MSE_independent)) ', p = ' num2str(P)]);
    
end


%% plot particular adjacency matrix

rep = 3;

taskNet = batch_log(1,1).taskNet(rep);
%tasksToPerform = batch_log(1,1).tasksToPerform(rep,:);
tasksToPerform = batch_log(1,1).tasksToPerform;
R_hidden = squeeze(batch_log(1,1).R_hidden(rep,:,:));
R_output = squeeze(batch_log(1,1).R_output(rep,:,:));


[pathwayCapacities maxCarryingCapacity  BK_MIS A_bipartite A_tasksIdx multiPerformance_mean multiPerformance_sem MSEdata A_dual A_dual_scores] = CogSci_validateMIS(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, multiCap);           

A_dual = A_dual_scores;

%%

% plot adjaciency matrix
figure(1);
map = colormap(copper);
map = map(fliplr(1:size(map,1)),:);
map(1,:) = 1;
colormap(map);

imagesc(A_dual); hold on;
hcb = colorbar;
%title('adjaciency matrix of interference graph', 'FontSize', 18);
set(gca, 'FontSize', 12);
set(gca, 'XTick', 1:size(A_dual,1));
set(gca, 'YTick', 1:size(A_dual,2));

for row = 1:size(A_dual,1)
    plot([row row]+0.5, [-1 size(A_dual,2)+1],'Color', [0.5 0.5 0.5]);
    plot([-1 size(A_dual,2)+1], [row row]+0.5,'Color', [0.5 0.5 0.5]);
end

colorTitleHandle = get(hcb,'Title');
titleString = 'score';
set(colorTitleHandle ,'String',titleString);
ylabel('tasks', 'FontSize', 18);
xlabel('tasks', 'FontSize', 18);

set(gcf, 'Color', [1 1 1])
set(gcf, 'Position', [100 100 400 300]);

%% PLOT TASK SIMILARITY MATRIX

rep = 3;

R_hidden = squeeze(batch_log(1,1).R_hidden(rep,:,:));
R_output = squeeze(batch_log(1,1).R_output(rep,:,:));

% threshold similarity matrix
R_hidden(R_hidden > corr_threshold) = 1;
R_output(R_output > corr_threshold) = 1;
R_hidden(R_hidden <= corr_threshold) = 0;
R_output(R_output <= corr_threshold) = 0;

ylimit = [-0.2 1];
A_dual = R_hidden;
% hidden similarity matrix
figure(1);
map = colormap(gray);
map = map(fliplr(1:size(map,1)),:);
map(1,:) = 1;
colormap(map);

imagesc(A_dual); hold on;
hcb = colorbar;
set(hcb, 'Ylim', ylimit);
set(gca, 'FontSize', 12);
set(gca, 'XTick', 1:size(A_dual,1));
set(gca, 'YTick', 1:size(A_dual,2));

for row = 1:size(A_dual,1)
    plot([row row]+0.5, [-1 size(A_dual,2)+1],'Color', [0.5 0.5 0.5]);
    plot([-1 size(A_dual,2)+1], [row row]+0.5,'Color', [0.5 0.5 0.5]);
end

colorTitleHandle = get(hcb,'Title');
titleString = 'r';
set(colorTitleHandle ,'String',titleString);
ylabel('tasks', 'FontSize', 18);
xlabel('tasks', 'FontSize', 18);
set(gcf, 'Color', [1 1 1])
set(gcf, 'Position', [100 100 310 250]);



A_dual = R_output;
% output similarity matrix
figure(2);
map = colormap(gray);
map = map(fliplr(1:size(map,1)),:);
map(1,:) = 1;
colormap(map);

imagesc(A_dual); hold on;
hcb = colorbar;
set(hcb, 'Ylim', ylimit);
set(gca, 'FontSize', 12);
set(gca, 'XTick', 1:size(A_dual,1));
set(gca, 'YTick', 1:size(A_dual,2));

for row = 1:size(A_dual,1)
    plot([row row]+0.5, [-1 size(A_dual,2)+1],'Color', [0.5 0.5 0.5]);
    plot([-1 size(A_dual,2)+1], [row row]+0.5,'Color', [0.5 0.5 0.5]);
end

colorTitleHandle = get(hcb,'Title');
titleString = 'r';
set(colorTitleHandle ,'String',titleString);
ylabel('tasks', 'FontSize', 18);
xlabel('tasks', 'FontSize', 18);
set(gcf, 'Color', [1 1 1])
set(gcf, 'Position', [100 100 310 250]);


%% PLOT MAXIMUM CAPACITY KINK

%% extract data 
rep = 20;

taskNet = batch_log(1,1).taskNet(rep);
performance = batch_log(1,1).multiPerformance_mean{20};
data1 = performance(2,:);
data2 = performance(2,:);

maxMIS = max(performance(1,:));

data1 = [];
data2 = [];

for cap = 2:NPathways % (maxMIS+1)
    
        input_MultiCap = multiCap{cap}.input;
        tasks_MultiCap = multiCap{cap}.tasks;
        train_MultiCap = multiCap{cap}.train;
        
        [~, ~, taskIdx] = unique(tasks_MultiCap, 'rows');
        
        [~, ~, MSE_multi] = taskNet.runSet(input_MultiCap, tasks_MultiCap, train_MultiCap);
            
        [GroupId,indx_i,index_j]=unique(taskIdx);
        GroupMean=arrayfun(@(k) mean(MSE_multi(index_j==k)),1:length(GroupId));
        
        data1 = [data1 min(GroupMean)];
        data2 = [data2 mean(MSE_multi)];
end

% plot

x = 2:NPathways;
plot(x, data2);

%% PLOT T-STATS DEPENDENT VS. INDEPENDENT

cap = 2;

good_MSE_data = [];
bad_MSE_data = [];

for rep = 1:length(batch_log(1,1).multiPerformance_mean)
   
    performance = batch_log(1,1).multiPerformance_mean{rep};
    if(any(performance(1,:) == cap))
       
        good_MSE_data = [good_MSE_data performance(2,performance(1,:)==cap)];
        bad_MSE_data = [bad_MSE_data performance(3,performance(1,:)==cap)];
        
    else 
       disp(['No cap=2 found at repetition = ' num2str(rep)]); 
    end
    
end

% plot

labelSize = 22;
lineWidth = 5;
fontName = 'Arial';

color = [50 50 50;150 150 150]/255; % nice color: 0 51 102
color = [0 0 100; 100 0 0]/255;

multiTest_good = mean(good_MSE_data);
multiTest_bad = mean(bad_MSE_data);

task_good_gen_sem = std(good_MSE_data)/sqrt(length(good_MSE_data));
task_bad_gen_sem = std(bad_MSE_data)/sqrt(length(bad_MSE_data));


bardata_mean = [mean(multiTest_good) mean(multiTest_bad)]/taskNet.Noutput;
            
bardata_sem = [task_good_gen_sem task_bad_gen_sem]/taskNet.Noutput;

%[bar_xtick ,hb,he] = errorbar_groups(bardata_mean', bardata_sem','bar_colors',color,'bar_width',2,'errorbar_width',0.3,'optional_errorbar_arguments',{'LineStyle','none','Marker','none','LineWidth',1.5});
hold on;
for i = 1:2
    bar(i, bardata_mean(i), 'FaceColor',color(i,:)); hold on;
    ha = errorbar(i, bardata_mean(i), bardata_sem(i), 'LineStyle','none', 'LineWidth', 0.25, 'Color', 'k'); 
end
%ylim([0 0.006]);
hb = get(ha,'children');
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
xleft = temp; xright = temp+1; 
Xdata(xleft) = Xdata(xleft) - .2;
Xdata(xright) = Xdata(xright) + .2;
set(hb(2),'Xdata',Xdata)

set(gca,'XTickLabel',{' ', ' '},'FontSize', labelSize, 'FontName', fontName,'FontWeight','normal');
ylabel('MSE','FontSize', labelSize, 'FontName', fontName,'FontWeight','normal');
%xlabel('tasks','FontSize', labelSize);
%title({'multitasking'},'FontSize', labelSize, 'FontName', fontName);
%ylim([0 0.07]);
xlim([0.5 2.5]);
set(gcf, 'Position', [800 800 400 250])
set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'independent', 'dependent'}) 
capacities = [batch_log(1, 1).MSEdata{1}.cap];

[H,P,CI,STATS] = ttest(bad_MSE_data - good_MSE_data,0,0.05,'right');
P

%% PLOT BIPARTITE GRAPH

%% plot adjacency matrix
rep = 3;

% find MIS
A_bipartite = batch_log(1,1).A_bipartite{rep};
%A_bipartite = A;

            W = 100;
            H = 100;
            Nlayers = 2;
            maxWeight = 1;
            maxLineWidth = 3; 
            
            Nhidden = size(A_bipartite,1);
            Noutput = size(A_bipartite,2);
            
            
 % init figure
            %clf;
            xlim([1 W]);
            ylim([1 H]);
            rectangle('Position',[0,0,W,H], 'FaceColor',[1 1 1], 'EdgeColor', [1 1 1])
            hold on;
            
            % calculate y-coordinates for all layers
            y_layer = [1:1:Nlayers] * H/(Nlayers+1);
            y_layer(1) = y_layer(1)-20;
            y_layer(end) = y_layer(end)+20;
            
            % calculate x-coordinates of units
            layer(1).x = [1:1:Nhidden] * W/(Nhidden+1);
            layer(2).x = [1:1:Noutput] * W/(Noutput+1);
            
            unitSize = [W/(Nhidden+1) W/(Noutput+1)] * 0.4;
            
                        % fill layer structure with y-coordinates
            for i = 1:length(layer)
               layer(i).y = repmat(y_layer(i),1,length(layer(i).x));
            end
            
            
            % draw connections between hidden & output units
            colors = distinguishable_colors(sum(sum(A_bipartite))+1);
            colors = 1-colors;
            colorIdx = 1;
            
            for hiddenUnit = 1:size(A_bipartite,1)
               for outputUnit = 1:size(A_bipartite,2)
                   lwidth = A_bipartite(hiddenUnit, outputUnit)/maxWeight * maxLineWidth;
                   draw = 1;
                   color = colors(1,:);
                   switch sign(lwidth)
                       case -1
                          color = '--k';
                       case 1
                           color = '-k';
                       case 0
                           draw = 0;
                   end
                   lwidth = abs(lwidth);
                   
                   color = colors(colorIdx,:);
                   if(lwidth > 0)
                      colorIdx = colorIdx + 1; 
                   end
                   color = [0 0 0];
                   
                   if(draw)
                        plot([layer(1).x(hiddenUnit) layer(2).x(outputUnit)], [layer(1).y(hiddenUnit) layer(2).y(outputUnit)], 'Color', color, 'LineWidth',lwidth*maxLineWidth-3);
                   end
               end
            end
            
            
            % draw units
            for layerIdx = 1:length(layer)
                for unit = 1:length(layer(layerIdx).x)
                    %plot(layer(layerIdx).x(unit), layer(layerIdx).y(unit), '.k', 'MarkerSize',unitSize)
                    %plot(layer(layerIdx).x(unit), layer(layerIdx).y(unit), '.w', 'MarkerSize',unitSize*0.95)
                    
                    x1 = layer(layerIdx).x(unit) - unitSize(layerIdx)/2;
                    y1 = layer(layerIdx).y(unit) - unitSize(layerIdx)/2;
                    x2 = unitSize(layerIdx);
                    y2 = unitSize(layerIdx);
                    
                    color = [1 1 1];
                    
                    color2 = [0 0 0];
                    lwidth = 1;
                    rectangle('Position',[x1,y1,x2,y2],'Curvature',[1,1], 'FaceColor',color, 'EdgeColor', color2, 'LineWidth', lwidth*maxLineWidth)
                end
            end
            
            set(gca, 'XTickLabels', {});
            set(gca, 'YTickLabels', {});
            set(gca,'Color',[1 1 1]);


%% EVALUATE DATA AT DIFFERENT CORRELATION THRESHOLD

corr_threshold = 0.85;
 
for rep = 1:length(batch_log(1,1).MSE_data)
    
            taskNet = batch_log(1,1).taskNet(rep);
            R_hidden = squeeze(batch_log(1,1).R_hidden(rep,:,:));
            R_output = squeeze(batch_log(1,1).R_output(rep,:,:));
            
            % validate MIS
            [pathwayCapacities maxCarryingCapacity  BK_MIS A_bipartite A_tasksIdx multiPerformance_mean multiPerformance_sem MSEdata A_dual A_dual_scores] = CogSci_validateMIS(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, multiCap);
            
            batch_log(1,1).maxCarryingCapacity(rep) = maxCarryingCapacity;
            batch_log(1,1).pathwayCapacities{rep} = pathwayCapacities;
            batch_log(1,1).BK_MIS{rep} = BK_MIS;
            batch_log(1,1).A_bipartite{rep} = A_bipartite;
            batch_log(1,1).A_dual{rep} = A_dual;
            batch_log(1,1).A_dual_scores{rep} = A_dual_scores;
            batch_log(1,1).A_tasksIdx{rep} = A_tasksIdx;
            batch_log(1,1).multiPerformance_mean{rep} = multiPerformance_mean;
            batch_log(1,1).multiPerformance_sem{rep} = multiPerformance_sem;
            batch_log(1,1).MSEdata{rep} = MSEdata;
    
end

%% SANITY CHECK: THRESHOLD PARAMETER - WEIGHT INITIALIZATION ROBUSTNESS

p_Matrix = NaN([size(batch_log) NPathways]);
capacities = 1:NPathways;


for init_scale_idx = 1:size(batch_log,1)
   
    for corr_thresh_idx = 1:size(batch_log, 2)
        
        MSE_data = nan(length(batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean), 2, NPathways);

        for rep = 1:length(batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean)

            for capIdx = 1:size(batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean{rep},2)

                MSE_data(rep, 1, batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean{rep}(1, capIdx)) = batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean{rep}(2, capIdx); % independent sets
                MSE_data(rep, 2, batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean{rep}(1, capIdx)) = batch_log(init_scale_idx,corr_thresh_idx).multiPerformance_mean{rep}(3, capIdx); % dependent sets

            end

        end
        
        MSE_dependent = zeros(1, size(MSE_data,1));
        MSE_independent = zeros(1, size(MSE_data,1));
        
        for cap = 2:NPathways
           
            MSE_independent = MSE_data(:,1,cap);
            MSE_dependent = MSE_data(:,2,cap);
            
            data = MSE_independent - MSE_dependent;
            data(isnan(data)) = [];
            
            [H,P,CI,STATS] = ttest(MSE_dependent - MSE_independent,0,0.05,'right');
            
            p_Matrix(init_scale_idx, corr_thresh_idx, cap) = P;
            
        end
            
    end
    
end

alpha = 0.05;

plot_Matrix = p_Matrix;
plot_Matrix(p_Matrix < alpha) = 1;
plot_Matrix(p_Matrix > alpha) = 0.5;
plot_Matrix(isnan(p_Matrix)) = 0;


%%

for cap = 2:NPathways
   
    plot_data = squeeze(plot_Matrix(:,:,cap));
    figure;
    imagesc(plot_data);
    map = [0 0 0; ...
           1 0 0; ...
           0 1 0];
    colormap(map);
    h = colorbar;
    set(h, 'YLim', [0 1]);
    title(['cap =  ' num2str(cap)], 'FontSize', 18);
    xlabel('green = OK, red = not OK, back = no data', 'FontSize', 18);
    set(gca, 'YTickLabel', init_scales);
    %set(gca, 'XTickLabel', corr_thresholds);
end

