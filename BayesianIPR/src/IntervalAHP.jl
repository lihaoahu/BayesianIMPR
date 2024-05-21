module IntervalAHP

using Turing
using LinearAlgebra
using Combinatorics
using ..BayesianIPR:IMPR,Interval

export iahp_analy,iahp,iahp_analy_uninfo,acceptability_index_mipr,uncertainty_index,iahp_posterior,rank_acceptability_index_mipr


"""
    iahp_analy(X,left,right,n,n_J,a,p,mu_prior,λ_mu,b_prior,λ_b)

根据区间偏好关系提供后验分布（解析解）

# Arguments
- `X`: 设计矩阵
- `left,right`: 左右端点向量
- `n,n_J`: 方案个数、偏好个数
"""
function iahp_analy(X,left,right,n,n_J,a,p,mu_prior,λ_mu,b_prior,λ_b)
    # 1. 计算对数排序向量后验均值
    m_mu = (λ_mu*I(n-1) + 2*X'*X)^(-1)*(X'*(left+right) + λ_mu*mu_prior);
    V_mu = (λ_mu*I(n-1) + 2*X'*X)^(-1);
    m_b = (ones(n_J)'*(right-left)+ λ_b*b_prior)/(λ_b+2*ones(n_J)'ones(n_J));
    v_b = 1/(λ_b+2*ones(n_J)'ones(n_J));

    alpha_phi = (2*n_J+p)/2;
    beta_phi = (a+left'*left+right'*right + λ_mu*mu_prior'*mu_prior + λ_b*b_prior^2  - m_mu'V_mu^(-1)*m_mu - m_b^2/v_b)/2;

    m_phi = alpha_phi/beta_phi;
    s_phi = alpha_phi/beta_phi^2;
    
    dof = 2*alpha_phi;

    S_mu = beta_phi/alpha_phi.*V_mu;
    C_mu = dof/(dof-2).*S_mu

    s_b = beta_phi/alpha_phi*v_b
    c_b = dof/(dof-2)*s_b

    M_hat = (2*X'*X)^(-1)*X'*(left+right)
    one = ones(Float64,n_J)
    b_hat = (2*one'*one)^(-1)*one'*(right-left)
    if λ_mu >0
        Ainv = ((2*X'*X)^(-1)+(λ_mu*I)^(-1))^(-1)
    else
        Ainv = λ_mu*I
    end
    
    if λ_b >0
        Binv = ((2*one'*one)^(-1)+(λ_b)^(-1))^(-1)
    else
        Binv = λ_b*I
    end

    return Dict(
        "ϕ" => Dict(
            "α" => alpha_phi,
            "β" => beta_phi,
            "mean" => m_phi,
            "var" => s_phi
        ),
        "μ" => Dict(
            "m" => m_mu,
            "V" => V_mu,
            "mean" => m_mu,
            "S" => S_mu,
            "Cov" => C_mu,
            "dof" => dof
        ),
        "b" => Dict(
            "m" => m_b,
            "v" => v_b,
            "mean" => m_b,
            "s" => s_b,
            "var" => c_b,
            "dof" => dof
        ),
        "pci" => beta_phi/(alpha_phi-1) * (n_J-1)/(n_J-n+1),
        "lci" => ((left-X*M_hat.+b_hat)'*(left-X*M_hat.+b_hat) + (right-X*M_hat.-b_hat)'*(right-X*M_hat.-b_hat))/(2*n_J+p-2) *(n_J-1)/(n_J-n+1),
        "prci" => (a+(mu_prior-M_hat)'*Ainv*(mu_prior-M_hat) + (b_prior-b_hat)^2*Binv)/(2*n_J+p-2) *(n_J-1)/(n_J-n+1)
    )
end

@model function iahp_posterior(post)
    ϕ ~ Gamma(post["ϕ"]["α"],1/post["ϕ"]["β"])
    σ² = 1/ϕ
    μ ~ MvNormal(post["μ"]["m"],σ²*round.(post["μ"]["V"];digits=6))

    b ~ Normal(post["b"]["m"],σ²*post["b"]["v"])
end

function iahp_analy_uninfo(X,left,right,n,n_J)

    return iahp_analy(X,left,right,n,n_J,0,0,zeros(n-1),0,0,0)
end

@model function iahp(X,left,right,n,n_J,a,p,u_prior,λ_mu,b_prior,λ_b)
    ϕ ~ Gamma(p/2,2/a)
    σ² = 1/ϕ
    μ ~ MvNormal(u_prior,σ²/λ_mu*I(n-1))

    bia ~ Normal(b_prior,sqrt(σ²/λ_b))

    left ~ MvNormal(X*μ .- bia,σ²*I(n_J))
    right ~ MvNormal(X*μ .+ bia,σ²*I(n_J))
end

function acceptability_index_mipr(chn::Chains)
    n = length(chn.name_map.parameters)-2+ 1
    chn = Array(chn);
    AI = zeros(Float32,(n,n))

    for t in 1:size(chn,1)
        u = cat(chn[t,2:n],[0],dims=1);
        v = exp.(u);
        w = v/sum(v);
        ow = sort(w,rev=true);
        for i in 1:n
            for k in 1:n
                if w[i] == ow[k]
                    AI[i,k] += 1
                end
            end
        end
    end
    return AI/size(chn,1)
end

function rank_acceptability_index_mipr(chn::Chains,P = [])
    n = length(chn.name_map.parameters)-2+ 1
    chn = Array(chn);
    if length(P) == 0
        P = permutations(1:n,n)
    end
    AI = Dict()
    for p in P
        AI[p] = 0
    end

    for t in 1:size(chn,1)
        u = cat(chn[t,2:n],[0],dims=1);
        v = exp.(u);
        w = v/sum(v);
        pw = sortperm(w,rev=true)
        for p in P
            p_bool = true
            for k in 2:length(p)
                if w[p[k]] > w[p[k-1]]
                    p_bool = false
                    break
                end
            end
            if p_bool
                AI[p] = AI[p] + 1
            end
        end
    end

    for p in P
        if AI[p] == 0
            delete!(AI,p)
        else
            AI[p] = AI[p]/size(chn,1)
        end
    end

    return AI
end

# function acceptability_index(R::MIPR,p::Real,a::Real,λ_mu::Real,λ_b::Real,rows::Vector{Int}=Vector{Int}([]),N::Int=10000)
#     Y = log.(R)
#     X,left,right,_,_,n,n_J = getdata(Y,rows);
#     M = iahp(X,left,right,n,n_J,a,p,λ_mu,λ_b)
#     chn = sample(M, NUTS(), N, progress=false)
#     return chn,acceptability_index(chn)
# end

function uncertainty_index(R::Matrix{Interval})
    s = 0
    for i in 1:size(R,1)
        for j in 1:size(R,1)
            s += log(R[i,j].up) - log(R[i,j].low)
        end
    end
    return s/size(R,1)/(size(R,1)-1)
end

include("IMPR/Tu_Wu.jl")

end