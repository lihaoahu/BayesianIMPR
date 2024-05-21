"""
    Interval(left,right)

返回左右端点分别为left、right的区间对象。

# Examples
```julia-repl
julia> A = Interval(0,1);

julia> A.low
```
"""
struct Interval
    low::Real
    up::Real
    Interval(low,up) = low > up ? error("非法区间！","[",low,",",up,"]") : new(low,up)
end

Base.log(A::Interval) = Interval(Base.log(A.low),Base.log(A.up))
Base.zero(Interval) = Interval(0,0)
Base.:+(A::Interval,b::Number) = Interval(A.low+b,A.up+b)
Base.:-(A::Interval,b::Number) = Interval(A.low-b,A.up-b)
Base.:/(A::Interval,d::Number) = Interval(min(A.low/d,A.up/d),max(A.low/d,A.up/d))
Base.:*(A::Interval,d::Number) = Interval(min(A.low*d,A.up*d),max(A.low*d,A.up*d))
Base.show(io::IO,A::Interval) = print(io,"[",round(A.low;digits=2),", ",round(A.up;digits=2),"]")
Base.one(Interval) = Interval(1,1)