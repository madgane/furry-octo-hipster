
user_indices = 10:10:100;

load savedData_plots;
A = reshape(plot_values,10,10);A(1,2) = A(1,3);

markerIndices = cell(5,1);

markerIndices{1,1} = 'b*-';
markerIndices{2,1} = 'gs-';
markerIndices{3,1} = 'rd-';
markerIndices{4,1} = 'co-';
markerIndices{5,1} = 'mp-';

hold on;
for i = 1:5
    plot(user_indices,A(:,i + 5),markerIndices{i,1});
end
hold off;
