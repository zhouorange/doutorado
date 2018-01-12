function [label_vector] = getLabelVector(label_vector,g_111,trainIter,N)

line_idx = N*(trainIter-1);
for y_col_idx=1:1:N
    column_idx = 0;
    line_idx = line_idx + 1;
    for y_line_idx=1:1:size(g_111,1)
        column_idx = column_idx + 1;
        label_vector(line_idx,column_idx) = real(g_111(y_line_idx,1));
        column_idx = column_idx + 1;
        label_vector(line_idx,column_idx) = imag(g_111(y_line_idx,1));
    end
end
