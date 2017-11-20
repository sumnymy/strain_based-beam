function variable_col=matrix2col(variable)
%To convert the variables from matrix to vector form
%----------------------------------------------------
%----------------------------------------------------
%Purpose:
%      To convert the variables from matrix to vector form
%Synopsis:
%      variable_col=matrix2col(variable)
%Variables definition:
%Input:
%      variable  -  variables in matrix form
% ----> h(nnode,12),pos(nnode,3),baseVector(nnode,9),strain(nelem,4)
%Output:
%      variable_col - variables in column vector form
% ---->h(12*nnode,1),pos(3*nnode,1),baseVector(9*nnode,1),strain(4*nelem,1)
%-------------------------------------------------------------------------
%get dimension of the variable
[r,c] = size(variable);
variable_col = zeros(c*r,1);
%convert matrix to vector
for i = 1:r
    for j = 1:c
        node = c*(i-1)+j;
        variable_col(node) = variable(i,j);
    end
end
