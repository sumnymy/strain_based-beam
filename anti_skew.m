%%written in 2017/9/24 by MY
%anti_skew---Convert a dual matrixx to its corresponding vector
function m=anti_skew(r)
m = zeros(3,1);
m(1)=r(3,2);
m(2)=r(1,3);
m(3)=r(2,1);

