%2017.9.20  by MY
%���ڽ�һ�������ת����һ�����������������һ������p2ת����һ������q2
%�ڴ˻����ϣ�����q1��q3
%R_mid��q Ϊת�������������pΪת��ǰ����������method=1��ʾ����ʵ�����method=2��ʾ��Сת������
%CF. D:\Research\2017.03---Ӧ�������Ԫ\doc\configuration.docx
function [q1,q2,q3]=R_mid(p1,p2,p3,q1,method)
%ȷ��Ϊ������
[r1,c1]=size(p1);
if r1<c1
    p1=p1';
end
[r2,c2]=size(p2);
if r2<c2
    p2=p2';
end
[r3,c3]=size(p3);
if r3<c3
    p3=p3';
end
[r4,c4]=size(q1);
if r4<c4
    q1=q1';
end
%%%p1---->q1
R = p1'*q1*eye(3)+skew(skew(p1)*q1)+skew(p1)*q1*(skew(p1)*q1)'/(1+p1'*q1);
q1=R*p1;
if method ==1
  q2=p2-p2'*q1*(p1+q1)/(1+p1'*q1);
  q3=p3-p3'*q1*(p1+q1)/(1+p1'*q1);
end

if method ==2
  q2=p2-p2'*q1*(p1+q1)/2;
  q3=p3-p3'*q1*(p1+q1)/2; 
end
