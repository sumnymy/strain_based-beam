function variable=col2matrix(variable_col,dimension)
%To convert the variables from matrix to vector form
%----------------------------------------------------
%----------------------------------------------------
%Purpose:
%      To convert the variables from vector to matrix form
%Synopsis:
%      variable=col2matrix(variable_col)
%Variables definition:
%Input:
%      variable_col - variables in column vector form
% ---->h(12*nnode,1),pos(3*nnode,1),baseVector(9*nnode,1),strain(4*nelem,1)
%      dimension   -  dimension of the variable at each node
% ---->h:12,pos:3,baseVector:9,strain:4
%Output:
%      variable  -  variables in matrix form
% ----> h(nnode,12),pos(nnode,3),baseVector(nnode,9),strain(nelem,4)
%-------------------------------------------------------------------------
%Initializing variable
len = length(variable_col);
r = len/dimension;
c = dimension;
variable = zeros(r,c);
%convert vector to matrix
for i = 1:r
    for j = 1:c
        node = c*(i-1)+j;
        variable(i,j) = variable_col(node);
    end
end