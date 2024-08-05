function getbit(n::Unsigned, i::Integer)
    return (n & (one(n) << i)) != 0
end

function setbit(n::Unsigned, i::Integer)
    return n | (one(n) << i)
end

function unsetbit(n::Unsigned, i::Integer)
    return n & ~(one(n) << i)
end

function flipbit(n::Unsigned, i::Integer)
    return n ⊻ (one(n) << i)
end

function flipbits(n::Unsigned, i::Integer, j::Integer)
    n = flipbit(n, i)
    n = flipbit(n, j)
    return n
end

function getspin(bool) :: Float64
    return -0.5 + bool
end

function getspin(n::Unsigned, i::Integer) :: Float64
    return getspin(getbit(n, i))
end