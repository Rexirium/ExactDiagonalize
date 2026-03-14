function readbit(bits::Int, pos::Int)::Bool
    return (bits >> (pos - 1)) & 1 == 1
end
# count number of ones between bit i and j (excluding i and j), i<j 
# and return the sign according to the odd or even of the number 
function signbetween(bits::Int, i::Int, j::Int)
    mask = 1<<(j-i-1) -1
    segbits = (bits>>i) & mask
    return (-1)^count_ones(segbits)
end

function flip(bits::Int, pos::Int)::Int
    return bits ⊻ (1 << (pos -1))
end

function splitbasis(bits::Int, b::Int)
    b >=0 || return 0, bits
    left = bits >> b
    right = bits & ((1<<b) - 1)
    return right, left
end

function numbitbasis(len::Int, num::Int)
    """
    generating all the L-bit integers with N bits read 1.
    """
    num > len && error("N is larger than L")
    num == 0 && return Int[0]
    basis = Int[]
    sizehint!(basis, binomial(len, num))
    maxind = (1 << len) - 1
    ind = (1 << num) - 1
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