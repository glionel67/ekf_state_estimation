function M = enforceMatrixSymmetry(M)
%enforceMatrixSymmetry: enforce the symmetry of the given square matrix M

assert(size(M, 1) == size(M, 2), 'enforceMatrixSymmetry: matrix must be square!');

M = 0.5 * (M + M');

% for i = 1 : size(M, 1)
%     for j = 1 : size(M, 2)
%         mij = 0.5 * M(i, j) + M(j, i);
%         M(i, j) = mij;
%         M(j, i) = mij;
%     end
% end

end
