export Tu_interval_priorites,Tu_consistency_index

function Tu_interval_priorites(R)
    B = log.(R.R);
    n = size(B,1);
    bl = map(x -> x.low,B);
    bu = map(x -> x.up,B);
    mu_w = mean(bl,dims=2) + mean(bu - bl,dims=2)/2 .- log(
        sum(
            exp.(
                mean(bl,dims=2) + mean(bu - bl,dims=2)/2
            )
        )
    )
    sigma_w = sqrt.(sum((bu - bl).^2,dims=2)/36/size(B,1)^2)
    return [Interval(exp(mu_w[i]-3*sigma_w[i]),exp(mu_w[i]+3*sigma_w[i])) for i in 1:n]
end

function Tu_consistency_index(R)
    B = log.(R.R)
    n = size(B,1);
    bl = map(x -> x.low,B);
    bu = map(x -> x.up,B);
    sigma2(i,j) = (
        n^2*(bu[i,j]-bl[i,j])^2
        + sum((bu[i,k] - bl[i,k])^2 for k in 1:n)
        + sum((bu[k,j] - bl[k,j])^2 for k in 1:n)
        )/36/n^2
    mu(i,j) = bl[i,j] - mean(bl[i,:]) + mean(bl[j,:]) + (bu[i,j]-bl[i,j])/2 - mean(bu[i,:]-bl[i,:])/2 + mean(bu[j,:]-bl[j,:])/2
    egci = 0
    for i in 1:n
        for j in 1:n
            if i<j
                egci += sigma2(i,j)+mu(i,j)^2
            end
        end
    end
    egci = 2*egci/(n-1)/(n-2)
    return egci
end
