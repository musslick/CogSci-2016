function MultiTasking_Result=Compute_MultiTask(Adjacency_BiPartite)

%
%
m = nnz(Adjacency_BiPartite);
[x y] = find(Adjacency_BiPartite);
%
A_dual = zeros(m);
%
for i=1:m
    for j = i+1:m
        first_in   = x(i);
        first_out  = y(i);
        second_in  = x(j);
        second_out = y(j);
        %
        if (first_in==second_in)
            A_dual(i,j) = 1;
        elseif (first_out==second_out)
            A_dual(i,j) = 1;
        elseif ~(isempty(find(y(find(x==first_in))==second_out))&&isempty(find(y(find(x==second_in))==first_out)))
            A_dual(i,j) = 1;
        end
    end
end
%
A_dual = (A_dual+A_dual');
%
MultiTask=findMIS(logical(A_dual),[1:m]);
%`
MultiTasking_Result = [x y MultiTask];