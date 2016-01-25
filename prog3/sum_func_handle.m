function sum = sum_func_handle(Xcol, Xcenter, fh_cell)
% Assistance function to realize summation over function handles for a function
% cell

dim = length(fh_cell);
n_col = size(Xcol,1);
n_cen = size(Xcenter,1);
sum = zeros(n_col,n_cen);
for i=1:dim
    sum = sum + fh_cell{i}(Xcol, Xcenter);
end

end