%2017.9.5  by MY
%���ڼ������㵥Ԫ�µ����ʣ��ο����ף�[1999][Crisfield]Objectivity of strain measures in the geometrically exact three-dimensional beam theory and its finite-element implementation.
%�������theta1,theta2�ֱ��ʾ��Ԫ���˵�����ת��������L��ʾ��Ԫ����
function k=curvature(theta1,theta2,L)
ta=0.5*(theta1+theta2);
td=theta2-theta1;
ta_sum=sqrt(sum(ta.^2));%ta��ģ��

k=1/L*(1/ta_sum^2*(1-sin(ta_sum)/ta_sum)*ta*ta'+sin(ta_sum)/ta_sum*eye(3)-(1-cos(ta_sum))/ta_sum^2*skew(ta))*td; %#ok<MHERM>