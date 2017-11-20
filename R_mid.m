%2017.9.20  by MY
%用于将一个坐标基转到另一个坐标基，满足其中一个向量p2转到另一个向量q2
%在此基础上，生成q1和q3
%R_mid，q 为转换后的向量基，p为转换前的向量基，method=1表示按真实情况，method=2表示按小转动近似
%CF. D:\Research\2017.03---应变基有限元\doc\configuration.docx
function [q1,q2,q3]=R_mid(p1,p2,p3,q1,method)
%确保为列向量
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
