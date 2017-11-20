%2017.9.5  by MY
%用于生成向量转动矩阵R， CF. D:\Research\2017.03---应变基有限元\doc\configuration.docx
%输入参数theta=(theta1;theta2;theta3)表示转轴向量；
%method用于选择矩阵生成方法：
%method=1: 式（2.13）
%method=2: 式（2.20）
%method=3：式（2.23）
%method=4:四元数，phi=2sin(theta/2)e

function R=rotation(theta,method)
%转动矩阵的基本形式
t_tidle = skew(theta);
t_sum = sqrt(sum(theta.^2));
%确保为列向量
[r,c]=size(theta);
if r<c
    theta=theta';
end
%Method1，原始方法
if method == 1
    R = eye(3)+sin(t_sum)/t_sum*t_tidle+(1-cos(t_sum))/sum(theta.^2)*t_tidle*t_tidle;
end
%FORM2，简化形式
if method == 2
    R = eye(3)+t_tidle+0.5*t_tidle*t_tidle;
end
%FORM3，简化形式
if method == 3
    R = eye(3)+(t_tidle+0.5*t_tidle*t_tidle)/(1+0.25*sum(theta.^2));
end
%四元数法，结果与方法1相同
if method == 4
    q0 = cos(t_sum/2);
    q = sin(t_sum/2)*theta/t_sum;
    R =(q0^2-q'*q)*eye(3)+2*q*q'+2*q0*skew(q);
end