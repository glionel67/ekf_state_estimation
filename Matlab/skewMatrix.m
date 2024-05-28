function mat = skewMatrix(vec)
%skewMatrix: convert the vector <vec> into its skew symmetric matrix form
% Args:
% - vec (vector [3x1]): vector
% Return:
% - mat (matrix [3x3]): skew symmetric matrix

mat = [0.0, -vec(3), vec(2);
    vec(3), 0.0, -vec(1);
    -vec(2), vec(1), 0.0];

end
