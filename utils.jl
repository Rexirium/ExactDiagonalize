function readbit(num::Int, pos::Int)
    return (num >> (pos -1)) & 1 == 1
end
# count number of ones between bit i and j (excluding i and j), i<j 
# and return the sign according to the odd or even of the number 
function signbetween(num::Int, i::Int, j::Int)
    mask = 1<<(j-i-1) -1
    segnum = (num>>i) & mask
    return (-1)^count_ones(segnum)
end

function splitbasis(num::Int, b::Int)
    b >=0 || return 0, num
    left = num >> b
    right = num & ((1<<b) - 1)
    return right, left
end

function numbitbasis(L::Int, N::Int)
    """
    generating all the L-bit integers with N bits read 1.
    """
    N > L && error("N is larger than L")
    N == 0 && return Int[0]
    basis = Int[]
    sizehint!(basis, binomial(L, N))
    maxind = (1<<L) - 1
    ind = (1<<N) - 1
    while ind <= maxind
        push!(basis, ind)
        u = ind & (-ind)
        v = ind + u
        next = v + ((v ⊻ ind) ÷ u) >> 2
        next > maxind && break
        ind = next
    end
    return basis
end