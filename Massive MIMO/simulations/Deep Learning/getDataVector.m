function [data_vector] = getDataVector(data_vector,Y1,trainIter,N)

line_idx = N*(trainIter-1);
for y_col_idx=1:1:size(Y1,2)
    column_idx = 0;
    line_idx = line_idx + 1;
    for y_line_idx=1:1:size(Y1,1)
        column_idx = column_idx + 1;
        data_vector(line_idx,column_idx) = real(Y1(y_line_idx,y_col_idx));
        column_idx = column_idx + 1;
        data_vector(line_idx,column_idx) = imag(Y1(y_line_idx,y_col_idx));
    end
end
