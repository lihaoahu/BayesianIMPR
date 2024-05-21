module AHP

using Turing
using LinearAlgebra

export ahp_analy,ahp
"""
ahp_analy(X,y,n,n_J,a,p,λ)

返回贝叶斯方法求解AHP问题的解析解

"""
function ahp_analy(X,y,n,n_J,a,p,λ)
    m_mu = (λ*I(n-1) + X'*X)^(-1)*(X'*y);
    S_mu = (λ*I(n-1) + X'*X)^(-1);

    alpha_phi = (n_J+p)/2;
    beta_phi = (a+y'*y-m_mu'S_mu^(-1)*m_mu)/2;

    m_phi = alpha_phi/beta_phi;
    s_phi = alpha_phi/beta_phi^2;

    e = (a+y'*y-m_mu'S_mu^(-1)*m_mu)/(n_J+p)
    m = [m_phi; m_mu]
    # v = [s_phi; e.*diag(S_mu)]

    return m,s_phi,e.*S_mu
end

@model function ahp(X,y,n,n_J,a,p,λ)
    ϕ ~ Gamma(p/2,2/a)
    σ² = 1/ϕ
    μ ~ MvNormal(zeros(Float64,n-1),σ²/λ*I(n-1))

    y ~ MvNormal(X*μ,σ²*I(n_J))
end
end