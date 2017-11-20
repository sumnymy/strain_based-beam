%%written in 2017/7/16 by MY
%skew---Convert a column vector to its dual matrix
function m_tilde=skew(m)
m_tilde = zeros(3,3);
m_tilde(1,2)=-m(3);
m_tilde(1,3)=m(2);
m_tilde(2,1)=m(3);
m_tilde(2,3)=-m(1);
m_tilde(3,1)=-m(2);
m_tilde(3,2)=m(1);
