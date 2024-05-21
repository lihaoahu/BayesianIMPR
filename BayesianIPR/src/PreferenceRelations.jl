abstract type PreferenceRelation{T} end

function Base.getindex(A::PreferenceRelation,I...)
    return Base.getindex(A.R,I...)
end

function Base.show(io::IO,::MIME"text/plain",A::PreferenceRelation{Interval})
    n = size(A.R,1)
    for i in 1:n
        for j in 1:n
            print(io,A.R[i,j],"\t")
        end
        println(io,"")
    end
end

struct MultiplicativePreferenceRelation{T} <: PreferenceRelation{T}
    R::Matrix{T}
end

struct AdditivePreferenceRelation{T} <: PreferenceRelation{T}
    R::Matrix{T}
end

NPR = NumericalPreferenceRelation = PreferenceRelation{Float64}
MPR = NumericalMultiplicativePreferenceRelation = MultiplicativePreferenceRelation{Float64}
IPR = IntervalPreferenceRelation = PreferenceRelation{Interval}
IMPR = MultiplicativeIntervalPreferenceRelation = MultiplicativePreferenceRelation{Interval}
IAPR = AdditiveIntervalPreferenceRelation = AdditivePreferenceRelation{Interval}


"""
    getdata(R,rows)

返回区间偏好关系矩阵上三角中指定行数据

# Arguments
- `R::Matirx{Interval}`: 区间性偏好关系矩阵
- `rows::Vector{Int}`: 所需行编号组成的向量
"""
function getdata(R::Matrix{Interval},rows::Vector{Int}=Vector{Int}([]))
    R = R
    n = size(R,1)
    if (length(rows)==0)
        n_J = sum(n-i for i in 1:n)
    else
        n_J = sum(n-i for i in 1:length(rows))
    end
    X = zeros(Int,(n_J,n-1));J = Vector{Vector{Int}}([]);l = zeros(Real,n_J);u = zeros(Real,n_J);r = zeros(Interval,n_J)
    k = 1
    for i in 1:n
        if length(rows)==0 || i in rows
            for j in 1:n
                x = zeros(Int,n-1)
                if i < n
                    x[i] = 1
                end
                if j < n
                    x[j] = -1
                end
                if !([i,j] in J) && !([j,i] in J) && i!=j
                    X[k,:] = x'
                    r[k] = R[i,j]
                    l[k] = R[i,j].low
                    u[k] = R[i,j].up
                    J = [J;[[i,j]]]
                    k += 1
                end
            end
        end
    end
    index = map(x->(!isnan(x.low)),r)
    X = X[index,:]
    l = l[index]
    u = u[index]
    J = J[index,:]
    r = r[index]
    n_J = Int(sum(index))
    # for i in 1:n_J
    #     if isnan.(r[i].low) || isnan(r[i].up)
    #         dele
    #     end
    # end

    return X,l,u,J,r,n,n_J
end

"""
    getall(R,rows)

返回区间偏好关系矩阵对角线外指定行数据

# Arguments
- `R::Matirx{Interval}`: 区间性偏好关系矩阵
- `rows::Vector{Int}`: 所需行编号组成的向量
"""
function getall(R::Matrix{Interval},rows::Vector{Int}=Vector{Int}([]))
    n = size(R,1)
    if (length(rows)==0)
        n_J = n*(n-1)
    else
        n_J = length(rows)*(n-1)
    end
    X = zeros(Int,(n_J,n-1));J = Vector{Vector{Int}}([]);l = zeros(Real,n_J);u = zeros(Real,n_J);r = zeros(Interval,n_J)
    k = 1
    for i in 1:n
        if length(rows)==0 || i in rows
            for j in 1:n
                x = zeros(Int,n-1)
                if i < n
                    x[i] = 1
                end
                if j < n
                    x[j] = -1
                end
                if i!=j
                    X[k,:] = x'
                    r[k] = R[i,j]
                    l[k] = R[i,j].low
                    u[k] = R[i,j].up
                    J = [J;[[i,j]]]
                    k += 1
                end
            end
        end
    end
    return X,l,u,J,r,n,n_J
end

function getdata(R::MPR)
    R = log.(R.R)
    n = size(R,1)
    nJ = sum(n-i for i in 1:n)
    X = zeros(Float64,(nJ,n-1));J = Vector{Vector{Int}}([]);y = zeros(Float64,nJ)
    k = 1
    for i in 1:n
        for j in 1:n
            x = zeros(Int,n-1)
            if i < n
                x[i] = 1
            end
            if j < n
                x[j] = -1
            end
            if !([i,j] in J) && !([j,i] in J) && i!=j
                X[k,:] = x'
                y[k] = R[i,j]
                J = [J;[[i,j]]]
                k += 1
            end
        end
    end
    index = map(x->(!isnan(x)),y)
    X = X[index,:]
    y = y[index]
    J = J[index,:]
    nJ = Int(sum(index))
    return X,y,n,nJ,J
end

function getdata(Rs::Vector{MPR})
    R1 = popfirst!(Rs)
    X,y,n,nJ,J = getdata(R1);
    for Ri in Rs
        Xi,yi,_,nJi,Ji = getdata(Ri);
        X = vcat(X,Xi);
        y = vcat(y,yi);
        J = vcat(J,Ji)
        nJ = nJ + nJi
    end
    return X,y,n,nJ,J
end

function getdata_group(Rs::Vector{MPR})
    R1 = popfirst!(Rs)
    X,y,n,nJ = getdata(R1);
    p = 1
    P = p*ones(Int,size(y))
    for Ri in Rs
        p += 1
        Xi,yi,_,nJi = getdata(Ri);
        X = vcat(X,Xi);
        y = vcat(y,yi);
        P = vcat(P,p*ones(Int,size(yi)))
        nJ = nJ + nJi
    end
    return X,y,P,n,nJ
end

function reorder(R::Union{IMPR,IAPR},order::Vector{Int}=Vector{Int}([]))
    if length(order) > 0
        n = length(order)
        New = zeros(eltype(R.R),(n,n))
        for i in 1:n
            for j in 1:n
                New[i,j] = R[order[i],order[j]]
            end
        end
        return typeof(R)(New)
    else
        return R
    end 
end