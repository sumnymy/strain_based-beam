%% 2017/7/17 by my
%Newton������������Է�����
%������������ֵ������4��P92
n = 2;
e1 = [1;0];
e2 = [0;1];
x=[0;0];           %��ֵ
eps = 10^(-12);    %����ˮƽ
M = 100;           %����������
k = 1;             %��������

while (k<=M)       %����������M
   h = [0.2;0.1];  
   F = -[3*x(1)-cos(x(1))-sin(x(2));4*x(2)-sin(x(1))-cos(x(2))];   %�����ұ�������
   J = [(3*(x(1)+h(1))-cos(x(1)+h(1))-sin(x(2))-3*x(1)+cos(x(1))+sin(x(2)))/h(1),...
       (3*x(1)-cos(x(1))-sin(x(2)+h(2))-3*x(1)+cos(x(1))+sin(x(2)))/h(2);...
       (4*x(2)-sin(x(1)+h(1))-cos(x(2))-4*x(2)+sin(x(1))+cos(x(2)))/h(1),...
       (4*(x(2)+h(2))-sin(x(1))-cos(x(2)+h(2))-4*x(2)+sin(x(1))+cos(x(2)))/h(2)];  %�ſ˱Ⱦ�����̵���ʽ
   deltax=J\F;
   if (norm(deltax)/norm(x)<=eps)    %������㾫�ȣ�������ѭ��
       break;
   else
       k=k+1;
       x=x+deltax;
   end
end

%%