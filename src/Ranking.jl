function rank_vector(x)
    sp = sortperm(x)
    k = zeros(Int64,length(sp))
    for i = 1:length(sp)
        k[sp[i]] = i
    end
    return k
end

function rank_to_binary(rank::Vector{Int64})
    n_candidates = length(rank)
    nu = triu!(ones(Int64, n_candidates, n_candidates), 1)
    nu = nu[rank, rank]
    return nu
end

function schulze(ranks)
    n_candidates, n_voters = size(ranks)
    D = sum([rank_to_binary(ranks[:, i]) for i = 1:n_voters])
    P = zeros(Int64, size(D))
    P[D .> D'] .= D[D .> D']
    
    for i = 1:n_candidates
        for j = 1:n_candidates
            if i != j
                for k = 1:n_candidates
                    P[j, k] = max(P[j, k], min(P[j, i], P[i, k]))
                end
            end
        end
    end
    return (n_candidates .- sum(P - P' .> 0, dims=2)[:]) |> sortperm |> invperm
end


function CV_array_to_rank(CV_total)
    rankings = hcat(rank_vector.(CV_total)...)
    return schulze(rankings) 
end
