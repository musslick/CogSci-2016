% author: Sebastian Musslick

function [pathwayCapacities maxCarryingCapacity  BK_MIS A_bipartite A_tasksIdx multiPerformance_mean multiPerformance_sem MSEdata A_dual A_dual_scores] = validateMIS(taskNet, tasksToPerform, R_hidden, R_output, corr_threshold, multiCap)

    sampleBadMulti = 1;

    NPathways = sqrt(taskNet.Ntask);

    % find MIS
    [pathwayCapacities maxCarryingCapacity  BK_MIS A_bipartite A_tasksIdx A_dual] = getMaxCarryingCapacity(R_hidden, R_output, corr_threshold);
    
    % formatting outputs
    BK_MIS = extend_BK_MIS(BK_MIS);
    A_tasksIdx(A_tasksIdx > 0) = tasksToPerform(A_tasksIdx(A_tasksIdx > 0));
    A_tasksOrder = A_tasksIdx(:);
    A_tasksOrder(A_tasksOrder == 0) = [];
    pathwayCapacities = [A_tasksOrder pathwayCapacities];

    % extract all possible multitasking conditions
    good_multiTaskConditions = zeros(size(BK_MIS,2), NPathways*NPathways);
    good_maximumCapacity = sum(BK_MIS);

    for multiCase = 1:size(BK_MIS,2)

        multiTaskIdx = pathwayCapacities(BK_MIS(:,multiCase) == 1,1);

        good_multiTaskConditions(multiCase,multiTaskIdx) = 1;

    end
    

    % extract all capacities (larger than just single task)
    allCapacities = unique(good_maximumCapacity);
    allCapacities(allCapacities == 1) = [];

    % for each multitasking capacity, calculate performance for
    % multitaskable tasks vs. non-multitaskable tasks
    multiPerformance_mean = zeros(3, length(allCapacities));
    multiPerformance_sem = zeros(size(multiPerformance_mean));
    multiPerformance_mean(1,:) = allCapacities;
    multiPerformance_sem(1,:) = allCapacities;
    
    MSEdata = {};

    for capIdx = 1:length(allCapacities)

        cap = allCapacities(capIdx);
        good_multiTaskConditions_cap = good_multiTaskConditions(good_maximumCapacity == cap,:);

        input_MultiCap = multiCap{cap}.input;
        tasks_MultiCap = multiCap{cap}.tasks;
        train_MultiCap = multiCap{cap}.train;

        goodMultiIdx = find(ismember(tasks_MultiCap, good_multiTaskConditions_cap, 'rows'));
        badMultiIdx = find(~ismember(tasks_MultiCap, good_multiTaskConditions_cap, 'rows'));
        
        if(sampleBadMulti) 
           [badTasksIDs,~,badTasksIdx] = unique(tasks_MultiCap(badMultiIdx,:),'rows');
           badTasksIdx_used = randsample(length(unique(badTasksIdx)),size(good_multiTaskConditions_cap,1));
           
           badMultiIdx_samp = badMultiIdx(ismember(badTasksIdx, badTasksIdx_used),:); 
        end

        % test good multitask performance
        [~, ~, MSE_multi_good] = taskNet.runSet(input_MultiCap(goodMultiIdx,:), tasks_MultiCap(goodMultiIdx,:), train_MultiCap(goodMultiIdx,:));
        [~,~,multi_good_tasksIdx] = unique(tasks_MultiCap(goodMultiIdx,:),'rows');
        
        % store mean performance on all good multitasking conditions
        [GroupId,~,index_j]=unique(multi_good_tasksIdx);
        MSE_multi_good_GroupMean=arrayfun(@(k) mean(MSE_multi_good(index_j==k)),1:length(GroupId));
        MSEdata(capIdx).cap = cap;
        MSEdata(capIdx).goodMSEData = MSE_multi_good_GroupMean;

        multiPerformance_mean(2, capIdx) = mean(MSE_multi_good);
        multiPerformance_sem(2,capIdx) = std(MSE_multi_good)/sqrt(length(MSE_multi_good));

        % test bad multitask performance

        [~, ~, MSE_multi_bad] = taskNet.runSet(input_MultiCap(badMultiIdx_samp,:), tasks_MultiCap(badMultiIdx_samp,:), train_MultiCap(badMultiIdx_samp,:));
        [~,~,multi_bad_tasksIdx] = unique(tasks_MultiCap(badMultiIdx_samp,:),'rows');
        
        % store mean performance on all good multitasking conditions
        [GroupId,~,index_j]=unique(multi_bad_tasksIdx);
        MSE_multi_bad_GroupMean=arrayfun(@(k) mean(MSE_multi_bad(index_j==k)),1:length(GroupId));
        MSEdata(capIdx).badMSEData = MSE_multi_bad_GroupMean;

        multiPerformance_mean(2, capIdx) = mean(MSE_multi_good);
        multiPerformance_sem(2,capIdx) = std(MSE_multi_good)/sqrt(length(MSE_multi_good));

        multiPerformance_mean(3, capIdx) = mean(MSE_multi_bad);
        multiPerformance_sem(3,capIdx) = std(MSE_multi_bad)/sqrt(length(MSE_multi_bad));
                                                    
    end
    
    % score A_dual graph
    A_dual_scores = 1-A_dual;
    
    % individual task scores
    A_dual_scores(1:(size(A_dual_scores,1)+1):end) = sum(BK_MIS,2);
    
    % task touples
    for row = 1:size(A_dual_scores,1)
        
        % cut lower triangular part
        %A_dual_scores(row, 1:(row-1)) = 0;
        
        for col = (row+1):size(A_dual_scores,2)
           if(A_dual_scores(row,col) == 1)
               touple_score = 0;
               for taskComp = 1:size(BK_MIS,2)
                  if(BK_MIS(row,taskComp) == 1 && BK_MIS(col,taskComp) == 1)
                     touple_score = touple_score+1;
                  end
               end
               A_dual_scores(row,col) = touple_score;
               A_dual_scores(col,row) = touple_score;
           end
        end
    end
    

end