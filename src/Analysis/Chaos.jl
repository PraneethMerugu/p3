module Chaos

using ..System
using ..Solvers
using ChaosTools
using DynamicalSystemsBase
using ComplexityMeasures
using DelayEmbeddings


export compute_lyapunov, compute_entropy, compute_lempel_ziv

"""
    compute_lyapunov(p::DuffingParams; T=1000.0)

Computes the maximum Lyapunov exponent using DynamicalSystems.jl.
"""
function compute_lyapunov(p::DuffingParams; T=1000.0)
    # Define CoupledODEs system
    ds = CoupledODEs(duffing!, [0.1, 0.1], p)
    
    # Compute MLE
    λ = lyapunov(ds, T)
    return λ
end

"""
    compute_entropy(trajectory_data)

Computes the Shannon entropy of the symbolic sequence derived from the trajectory.
Partition: x < 0 (0), x >= 0 (1).
"""
function compute_entropy(x_data::Vector{Float64})
    # Symbolize
    symbols = [x < 0 ? 0 : 1 for x in x_data]
    
    # Compute probabilities
    counts = Dict{Int, Int}()
    for s in symbols
        counts[s] = get(counts, s, 0) + 1
    end
    
    total = length(symbols)
    probs = [c/total for c in values(counts)]
    
    # Shannon Entropy
    H = -sum(p * log2(p) for p in probs if p > 0)
    return H
end

"""
    compute_lempel_ziv(trajectory_data)

Computes the normalized Lempel-Ziv complexity of the symbolic sequence.
"""
function compute_lempel_ziv(x_data::Vector{Float64})
    # Symbolize
    symbols = [x < 0 ? 0 : 1 for x in x_data]
    
    # Manual LZ76 implementation to avoid dependency issues
    n = length(symbols)
    if n <= 1
        return 0.0
    end
    
    c = 1
    l = 1
    i = 0
    k = 1
    k_max = 1
    
    while c + k <= n
        if symbols[c+k] == symbols[l+k]
            k += 1
            if l + k > c
                c += k
                i += 1
                l = 1
                k = 1
                k_max = 1
            else
                k_max = max(k_max, k)
            end
        else
            l += 1
            k = 1
            if l == c
                c += k_max
                i += 1
                l = 1
                k_max = 1
            end
        end
    end
    
    C = i + 1
    
    # Normalization: C_norm = C / (N / log2(N))
    return C / (n / log2(n))
end

end # module Chaos
