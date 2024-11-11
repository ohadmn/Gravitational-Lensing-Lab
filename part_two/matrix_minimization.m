function[x, y] = matrix_minimization(mat)
    min = flintmax;
    column = 1;
    row = 1;
    mat_dimensions_vector = size(mat);
    mat_dimension = mat_dimensions_vector(1);
    flag = false;
    f_selector = 1;
    while flag == false
        flag=true;
        if f_selector == 1
            for i = 1 : mat_dimension
                if min > mat(row,i)
                    flag=false;
                    min = mat(row,i);
                    column = i;
                end
            end
            f_selector = 2;
        else
            flag=true;
            for j = 1 : mat_dimension
                if min > mat(j,column)
                    min = mat(j,column);
                    flag=false;
                    row = j;
                end
            end
            f_selector = 1;
        end
    end
    x = column;
    y = row;
end