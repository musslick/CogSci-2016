clear; clc; close all;
%
%
max_NetWork_Size = 50;
time_Profile = NaN(1,max_NetWork_Size);
%
MultiTask_MaX = NaN(max_NetWork_Size);
MultiTask_MiN = NaN(max_NetWork_Size);
%
%%
%
for current_N = 1:max_NetWork_Size
    %
    tic
    %
    parfor path_Overlap = 1:current_N
        %
        %
        irrelevant_Coonection = path_Overlap-1;
        %
        temp_A22 = eye(current_N - (irrelevant_Coonection));
        %
        temp_A21_factor = zeros(current_N - (irrelevant_Coonection),1);
        temp_A21_factor(1) = 1;
        temp_A21 = temp_A21_factor*ones(1,irrelevant_Coonection);
        %
        temp_A1 = ones(irrelevant_Coonection,current_N);
        %
        A_bipartite_minDC = [temp_A1; [temp_A21 temp_A22]]';
        MT_Cabability_MaX = Compute_MultiTask(A_bipartite_minDC);
        %
        MultiTask_MaX(current_N,path_Overlap) = sum(MT_Cabability_MaX(:,3));
        %
        %
        %
        %
        row_first = [ones(1,path_Overlap) zeros(1,current_N-path_Overlap)];
        col_first = circshift([zeros(1,current_N-path_Overlap) ones(1,path_Overlap)]',1);
        %
        A_bipartite_maxDC = toeplitz(col_first,row_first);
        MT_Cabability_MiN = Compute_MultiTask(A_bipartite_maxDC);
        %
        MultiTask_MiN(current_N,path_Overlap) = sum(MT_Cabability_MiN(:,3));
        %
        clc
    end
    %
    time_Profile(current_N) = toc;
    disp(['NetWork Size: ',num2str(current_N),' ::: Time lapsed - ',num2str(time_Profile(current_N))])
    %
end
%
%
%%
[X_Axis,Y_Axis] = meshgrid(1:max_NetWork_Size);
% X_Axis: pathway overlap (P)
% Y_Axis: network size (N)
choice_Complexity = (Y_Axis-X_Axis+ones(max_NetWork_Size)).*(log2(X_Axis) + log2(Y_Axis));
range = (MultiTask_MaX-MultiTask_MiN);
%
%
figure(1)
colormap('jet')
%
h1 = surf(X_Axis,Y_Axis,MultiTask_MaX,choice_Complexity,'EdgeColor','none','FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',MultiTask_MaX);
view(20,15)
%
ylabel('Network Size - N'); xlim([1 1+max_NetWork_Size]);
xlabel('Path Overlap - P'); ylim([1 1+max_NetWork_Size]);
zlabel('Multitasking Capability with Min. DC');
alphamap(linspace(0.5,1,64));
%
hold on
for i = 1:5
    plot3(Y_Axis(:,i*10),X_Axis(:,i*10),MultiTask_MaX(i*10,:),'b','linewidth',2);
end
%
%
figure(2)
colormap('jet')
%
h2 = surf(X_Axis,Y_Axis,MultiTask_MiN,choice_Complexity,'EdgeColor','none','FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',MultiTask_MaX);
view(20,15)
%
ylabel('Network Size - N'); xlim([1 1+max_NetWork_Size]);
xlabel('Path Overlap - P'); ylim([1 1+max_NetWork_Size]);
zlabel('Multitasking Capability with Max. DC');
alphamap(linspace(0.5,1,64));
%
hold on
for i = 1:5
    plot3(Y_Axis(:,i*10),X_Axis(:,i*10),MultiTask_MiN(i*10,:),'b','linewidth',2);
end
%
%
figure(3)
colormap('jet')
%
h3 = surf(X_Axis,Y_Axis,range,choice_Complexity,'EdgeColor','none','FaceAlpha','flat','AlphaDataMapping','scaled','AlphaData',MultiTask_MaX);
view(20,15)
%
ylabel('Network Size - N'); xlim([1 1+max_NetWork_Size]);
xlabel('Path Overlap - P'); ylim([1 1+max_NetWork_Size]);
zlabel('Range of Multitasking Capability');
alphamap(linspace(0.5,1,64));
%
hold on
for i = 1:5
    plot3(Y_Axis(:,i*10),X_Axis(:,i*10),range(i*10,:),'b','linewidth',2);
end
%
%
colormap(jet)
%
h2 = surf(X_Axis,Y_Axis,MultiTask_MiN);
shading interp
alpha(0.7)
view(20,15)
%
ylabel('Network Size - N'); xlim([1 1+max_NetWork_Size]);
xlabel('Path Overlap - P'); ylim([1 1+max_NetWork_Size]);
zlabel('Multitasking Capability with Max. DC');
%
hold on
for i = 1:5
    plot3(Y_Axis(:,i*10),X_Axis(:,i*10),MultiTask_MiN(i*10,:),'k','linewidth',2);
end
view(20,15)