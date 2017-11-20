%2017.9.5  by MY
%������������ת������R�� CF. D:\Research\2017.03---Ӧ�������Ԫ\doc\configuration.docx
%�������theta=(theta1;theta2;theta3)��ʾת��������
%method����ѡ��������ɷ�����
%method=1: ʽ��2.13��
%method=2: ʽ��2.20��
%method=3��ʽ��2.23��
%method=4:��Ԫ����phi=2sin(theta/2)e

function R=rotation(theta,method)
%ת������Ļ�����ʽ
t_tidle = skew(theta);
t_sum = sqrt(sum(theta.^2));
%ȷ��Ϊ������
[r,c]=size(theta);
if r<c
    theta=theta';
end
%Method1��ԭʼ����
if method == 1
    R = eye(3)+sin(t_sum)/t_sum*t_tidle+(1-cos(t_sum))/sum(theta.^2)*t_tidle*t_tidle;
end
%FORM2������ʽ
if method == 2
    R = eye(3)+t_tidle+0.5*t_tidle*t_tidle;
end
%FORM3������ʽ
if method == 3
    R = eye(3)+(t_tidle+0.5*t_tidle*t_tidle)/(1+0.25*sum(theta.^2));
end
%��Ԫ����������뷽��1��ͬ
if method == 4
    q0 = cos(t_sum/2);
    q = sin(t_sum/2)*theta/t_sum;
    R =(q0^2-q'*q)*eye(3)+2*q*q'+2*q0*skew(q);
end