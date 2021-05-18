function [m_inf,m_inf_prime] = m_infty(u,par)

m_inf = 0.5.*(1 + tanh((u - par.u1)./par.u2));


if nargout > 1
    
    m_inf_prime = 0.5.*(sech((u - par.u1)./par.u2).^2)./par.u2;
    
end
