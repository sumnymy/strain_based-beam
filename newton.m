%% 2017/7/17 by my
%Newton迭代法求非线性方程组
%算例来自于数值分析第4版P92
n = 2;
e1 = [1;0];
e2 = [0;1];
x=[0;0];           %初值
eps = 10^(-12);    %精度水平
M = 100;           %最大迭代次数
k = 1;             %迭代次数

while (k<=M)       %最大迭代次数M
   h = [0.2;0.1];  
   F = -[3*x(1)-cos(x(1))-sin(x(2));4*x(2)-sin(x(1))-cos(x(2))];   %方程右边列向量
   J = [(3*(x(1)+h(1))-cos(x(1)+h(1))-sin(x(2))-3*x(1)+cos(x(1))+sin(x(2)))/h(1),...
       (3*x(1)-cos(x(1))-sin(x(2)+h(2))-3*x(1)+cos(x(1))+sin(x(2)))/h(2);...
       (4*x(2)-sin(x(1)+h(1))-cos(x(2))-4*x(2)+sin(x(1))+cos(x(2)))/h(1),...
       (4*(x(2)+h(2))-sin(x(1))-cos(x(2)+h(2))-4*x(2)+sin(x(1))+cos(x(2)))/h(2)];  %雅克比矩阵差商的形式
   deltax=J\F;
   if (norm(deltax)/norm(x)<=eps)    %如果满足精度，则跳出循环
       break;
   else
       k=k+1;
       x=x+deltax;
   end
end

%%