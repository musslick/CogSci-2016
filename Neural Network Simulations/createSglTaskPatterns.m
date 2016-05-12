function [inputSgl tasksSgl trainSgl tasksIdxSgl] = createSglTaskPatterns(NPathways, Nactive, NFeatures)
% createTrainingPatterns generate single tasking and multitasking training
% patterns
% 
% PARAMETERS:
%
% NPathways     ...Number of pathways (i.e. number of feature dimensions % output dimensions)
% Nactive       ...How many tasks to perform simultaneously
% NFeatures     ...lists all feature dimensions
%
% RETURN VALUES:
%
% inputSgl      ...Single-task training set for input stimuli 
% tasksSgl      ...Single-task training set for input tasks
% trainSgl      ...Single-task training set for correct network outputs
%
% inputMulti    ...Multi-task training set for input stimuli 
% tasksMulti    ...Multi-task training set for input tasks
% trainMulti    ...Multi-task training set for correct network outputs
%

% specify number of possible task combinations (this number blows up as we
% increase NPathways) for multitasking condition (2 tasks at the same time)
% also specify the number of possible input stimuli (only one feature per
% feature dimension can be active)

% author: Sebastian Musslick

stimCombs = [1:NFeatures];
for i = 2:NPathways
    stimCombs = combvec(stimCombs,[1:NFeatures]); 
end
currFeaturesDims = 1:NPathways;
% create the task environment for performing just one task (single task)

% loop through each new training set
input_sglT = [];
tasks_sglT = [];
train_sglT = [];
tasksIdx_sglT = [];

% number of input patterns 
Nsets_sglT = NPathways*NPathways * size(stimCombs,2);

for currT = 1:(NPathways*NPathways)
    
    % build task input
    currTasks = zeros(1,NPathways*NPathways);
    currTasks(currT) = 1;
    currTasks = repmat(currTasks(:)', size(stimCombs,2), 1); % backtransform: reshape(currTasks(1,:),10,10)'

    currTasksM = reshape(currTasks(1,:),NPathways,NPathways)';
    
    % build feature input
    currInput = zeros(size(stimCombs,2),NPathways*NFeatures);
    
    for i = 1:size(currInput,1)
        currInput(i,(currFeaturesDims-1).*(NFeatures)+stimCombs(:,i)') = 1;
    end
    
    % build output
    currTrain = zeros(size(stimCombs,2),NPathways*NFeatures);
    
    [relevantInput,relevantOutput] = find(currTasksM ==1);
    for i = 1:size(currInput,1)
        currTrain(i,[((relevantOutput-1)*NFeatures+1):((relevantOutput-1)*NFeatures+NFeatures)]) = currInput(i,((relevantInput-1)*NFeatures+1):((relevantInput-1)*NFeatures+NFeatures));
    end

    % build full training set
    tasks_sglT = [tasks_sglT; currTasks];
    input_sglT = [input_sglT; currInput];
    train_sglT = [train_sglT; currTrain];
    
    curr_tasksIdx = currT*ones(size(stimCombs,2),1);
    tasksIdx_sglT = [tasksIdx_sglT; curr_tasksIdx];
    
end

inputSgl = input_sglT;
tasksSgl = tasks_sglT;
trainSgl = train_sglT;
tasksIdxSgl = tasksIdx_sglT;

end