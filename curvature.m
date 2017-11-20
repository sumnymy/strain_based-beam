%2017.9.5  by MY
%用于计算两点单元下的曲率，参考文献：[1999][Crisfield]Objectivity of strain measures in the geometrically exact three-dimensional beam theory and its finite-element implementation.
%输入参数theta1,theta2分别表示单元两端的总体转动向量；L表示单元长度
function k=curvature(theta1,theta2,L)
ta=0.5*(theta1+theta2);
td=theta2-theta1;
ta_sum=sqrt(sum(ta.^2));%ta的模长

k=1/L*(1/ta_sum^2*(1-sin(ta_sum)/ta_sum)*ta*ta'+sin(ta_sum)/ta_sum*eye(3)-(1-cos(ta_sum))/ta_sum^2*skew(ta))*td; %#ok<MHERM>