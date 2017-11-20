%2017.9.5 by MY
%用于通过转动矩阵R生成转轴向量v， CF. D:\Research\2017.03---应变基有限元\doc\configuration.docx
%说明:用不同方法生成的psedo_vector是不同的，但都是成比例的
%method=1：转动矩阵直接生成
%method=2：四元素表示法生成
function v=psedo_vector(R,method)
tic;
if method == 1
    Ra = 0.5*(R-R');
    sv = 4*Ra/(1+trace(R));
    v = [sv(3,2);sv(1,3);sv(2,1)];
end
if method == 2
    if trace(R)>=max([R(1,1),R(2,2),R(3,3)])
        a=trace(R);
        q0=0.5*sqrt(1+a);
        q = [R(3,2)-R(2,3);R(1,3)-R(3,1);R(2,1)-R(1,2)]/(4*q0);
        v = 2*q*acos(q0)/sqrt(1-q0^2);
    else
        q=zeros(3,1);
        [a,i]=max([R(1,1),R(2,2),R(3,3)]);
        q(i)=sqrt(0.5*a+0.25*(1-trace(R)));
        switch i
            case 1
            q0=0.25*(R(3,2)-R(2,3))/q(i);
            q(2)=0.25*(R(1,2)+R(2,1))/q(i);
            q(3)=0.25*(R(1,3)+R(3,1))/q(i);
            case 2
            q0=0.25*(R(1,3)-R(3,1))/q(i);
            q(1)=0.25*(R(2,1)+R(1,2))/q(i);
            q(3)=0.25*(R(3,2)+R(2,3))/q(i);
            case 3
            q0=0.25*(R(2,1)-R(1,2))/q(i);
            q(2)=0.25*(R(2,3)+R(3,2))/q(i);
            q(1)=0.25*(R(1,3)+R(3,1))/q(i);    
        end
        v = 2*q*acos(q0)/sqrt(1-q0^2);
    end
end
toc
