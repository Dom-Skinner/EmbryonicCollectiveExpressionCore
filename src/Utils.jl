using LinearAlgebra
using StatsBase
using Random
using MultivariateStats
using Distributed

function vech(A::AbstractMatrix{T}) where T
    # Steven Johnson code
    m = LinearAlgebra.checksquare(A)
    v = Vector{T}(undef,(m*(m+1))>>1)
    k = 0
    for j = 1:m, i = j:m
        @inbounds v[k += 1] = A[i,j]
    end
    return v
end

function avech(A::AbstractMatrix{T}) where T
    # Steven Johnson code
    m = LinearAlgebra.checksquare(A)
    v = Vector{T}(undef,(m*(m-1))>>1)
    k = 0
    for j = 1:m, i = (j+1):m
        @inbounds v[k += 1] = A[i,j]
    end
    return v
end

function ivech(v::AbstractVector{T}) where T
    m = Int64(round( sqrt(2*length(v)+0.25)-0.5))
    A = Matrix{T}(undef,m,m)
    k = 0
    @inbounds for j = 1:m, i = j:m
        k += 1
        # A[i,j] = v[k] # only assign upper part of matrix
         A[j,i] = v[k]
    end
    return Symmetric(A)
end

function iavech(v::AbstractVector{T}) where T
    m = Int64(round( sqrt(2*length(v)+0.25)+0.5))
    A = zeros(T,m,m)
    k = 0
    @inbounds for j = 1:m, i = (j+1):m
        k += 1
         A[i,j] =  v[k]
         A[j,i] = -v[k]
    end
    return A
end

function pos_def_cone(A::AbstractMatrix{T}) where T
    λ,e = eigen(A)
    return sum( [(λ[i] > 0 ? λ[i] : 0)  * (e[:,i] * e[:,i]') for i = 1:length(λ)])
end



function indicator(t,type_)
    ind = zeros(size(type_))
    ind[findall(type_.==unique(type_)[t])] .= 1
    return ind
end

entropy(Mod) = 0.5*size(Mod.d.Σ,1)*(log(2*pi)+1)+ 0.5*log(det(Mod.d.Σ))

function permute_mat(X,s)
    Xp = copy(X)
    idx = randperm(size(X,1))[1:s]
    for i in idx
        Xp[i,:] .= shuffle(Xp[i,:])
    end
    return Xp,idx
end

@inline function triangle_inv(i,j,n)
    return (i > j) ? Int64(1 + (i-j) + n*(n+1)/2 - (n-j+2)*(n-j+1)/2) : Int64(1 + (j-i) + n*(n+1)/2 - (n-i+2)*(n-i+1)/2)
end

function setdiff_with_repeats(x::Array{T},y::Array{T}) where T
    # This function is  a setdiff like-function, but respects the number of
    # each occurence. E.g. setdiff([1],[1,1,2]) = [2] but 
    # setdiff_with_repeats([1],[1,1,2]) = [1,2]. Implementation not performant.
    to_return = []
    xs = countmap(x)
    ys = countmap(y)
    
    total_keys = vcat(collect(keys(xs)),collect(keys(ys))) |> unique

    for k in total_keys
        if !haskey(xs,k)
            xs[k] = 0
        end
        if !haskey(ys,k)
            ys[k] = 0
        end
        push!(to_return,repeat([k],abs(xs[k] - ys[k])))
    end
    return vcat(to_return...)
end


function arrange_GE_data(PC_GE::PCA_GeneExpression;rand_sample=true,bootstrap=false)
    @assert !ismissing(PC_GE.GE.cell_type_labels)

    if bootstrap
        return arrange_GE_data(PCA_GeneExpression(embryo_bootstrap(PC_GE.GE),PC_GE.nPC,PC_GE.M,
            PC_GE.rem_cluster_mean,PC_GE.mother_mean,PC_GE.time_mean),
                    rand_sample=rand_sample,bootstrap=false)
    end
    nc  = PC_GE.GE.stage
    nPC = PC_GE.nPC
    embryos = unique(PC_GE.GE.embryo_id)
    local_clusters = [PC_GE.GE.cell_type_labels[PC_GE.GE.embryo_id .== e] for e in embryos]
    ref_cluster = get_type(PC_GE.GE)

    
    X,_ = remove_variance(PC_GE.GE,rem_cluster_mean=PC_GE.rem_cluster_mean,mother_mean=PC_GE.mother_mean,
        time_mean=PC_GE.time_mean,randomize=rand_sample)
    X = MultivariateStats.predict(PC_GE.M,X)

    Z = [fill(NaN, nPC, nc) for idx in 1:length(embryos)]

    for i in 1:length(embryos)
        cond = findall(PC_GE.GE.embryo_id .== embryos[i])
   
        Z[i][:,1:length(cond)] = X[:,cond]
   
        # Rearrange to get in cell type order
        local_cluster_tmp = vcat(local_clusters[i],setdiff_with_repeats(local_clusters[i],ref_cluster))
        Z[i] .= Z[i][:,sortperm(local_cluster_tmp)]        
    end
    
    return Z
end

function arrange_GE_data(GE::GeneExpressionStage;rand_sample=true,bootstrap=false)
    @assert !ismissing(GE.cell_type_labels)

    nc  = GE.stage
    ng = GE.ng
    embryos = unique(GE.embryo_id)
    local_clusters = [GE.cell_type_labels[GE.embryo_id .== e] for e in embryos]
    ref_cluster = get_type(GE)


    Z = [fill(NaN, ng, nc) for idx in 1:length(embryos)]
    count_sample = copy(GE.counts)
    

    if bootstrap
        idx_arr = rand(1:length(embryos),length(embryos))
    else
        idx_arr = 1:length(embryos)
    end
    

    for i in 1:length(embryos)
        cond = findall(GE.embryo_id .== embryos[idx_arr[i]])
        if rand_sample
            Z[i][:,1:length(cond)] = count_sample[:,cond] .+ GE.error_bars[:,cond].*randn(size(GE.counts,1),length(cond))
        else
            Z[i][:,1:length(cond)] = count_sample[:,cond]
        end
        # Rearrange to get in cell type order
        local_cluster_tmp = vcat(local_clusters[idx_arr[i]],setdiff_with_repeats(local_clusters[idx_arr[i]],ref_cluster))
        Z[i] .= Z[i][:,sortperm(local_cluster_tmp)]        
    end
    
    return Z
end


function find_ng(GE::GeneExpressionStage)
    return GE.ng
end
function find_ng(PC_GE::PCA_GeneExpression)
    return PC_GE.nPC
end


function find_keep_vars(missing_cells,ref_cluster,stage)
    rem = []
    for c in missing_cells
        choices = findall(isequal(c),ref_cluster)
        for x in choices
            if x ∈ rem
                continue
            else
                push!(rem,x)
                break
            end
        end
    end
    return setdiff(1:stage,rem)
end
function generate_synthetic_data(GE::GeneExpressionStage,Mod::MaxEntModel;K=1,nPCA=0)
    synth_data = rand_sample(Mod,K*GE.N_samples)

    if nPCA > 0
        @assert Mod.ng == nPCA
        X,subtracted = remove_variance(GE,rem_cluster_mean=true,mother_mean=false,time_mean=true,randomize=false)
        M = fit(PCA, X; maxoutdim=nPCA);
        subtracted = repeat(subtracted,1,K)
        synth_data = [reconstruct(M,s) for s in synth_data]
    end


    ref_cluster = get_type(GE)
    synth_data_array = zeros(size(GE.counts,1),K*size(GE.counts,2))
    synth_err_bars = repeat(GE.error_bars,1,K)
    embryo_id = vcat([GE.embryo_id .+ (i-1)*(1 + maximum(GE.embryo_id) - minimum(GE.embryo_id)) for i = 1:K]...)
    cell_type_labels = repeat(GE.cell_type_labels,K)
    
    for (i,e) in enumerate(unique(embryo_id))
        cond = findall(isequal(e),embryo_id)
        sp = sortperm(cell_type_labels[cond])
        cell_type_labels[cond] .= cell_type_labels[cond][sp]
        synth_err_bars[:,cond] .= synth_err_bars[:,cond][:,sp] 
        local_cluster = cell_type_labels[cond]
        
        missing_cells =setdiff_with_repeats(local_cluster,ref_cluster)
        keep_vars =find_keep_vars(missing_cells,ref_cluster,GE.stage)
        synth_data_array[:,cond] .= synth_data[i][:,keep_vars]
        if nPCA > 0
            synth_data_array[:,cond] .+= subtracted[:,cond][:,sp] 
        end
    end


    return GeneExpressionStage(GE.stage,"synthetic",GE.genename,embryo_id,
        synth_data_array,cell_type_labels,GE.cell_type_sizes,synth_err_bars,GE.ng,GE.N_samples)
end


function embryo_bootstrap(GE::GeneExpressionStage)
    embryos = unique(GE.embryo_id)
    idx_arr = rand(1:length(embryos),length(embryos))
    new_idx = vcat([findall(GE.embryo_id .== embryos[k]) for k in idx_arr]...)
    new_id = vcat([i*ones(sum(GE.embryo_id .== embryos[idx_arr[i]])) for i in 1:length(idx_arr)]...)
    return GeneExpressionStage(GE.stage,"bootstrap",GE.genename,new_id,
    GE.counts[:,new_idx],GE.cell_type_labels[new_idx],GE.cell_type_sizes,GE.error_bars[:,new_idx],GE.ng,GE.N_samples)
end




function remove_variance(GE;rem_cluster_mean=true,mother_mean=false,time_mean=true,randomize=false)
    # Function that removes sources of variance, whether cluster specific mean or temporal projection 
    if randomize
        X = (GE.counts.+ GE.error_bars.*randn(size(GE.error_bars)))
    else
        X = copy(GE.counts)
    end

    subtracted = zeros(size(X))

    if rem_cluster_mean
        c_mean = cluster_mean(X,GE.cell_type_labels)
        X .-= c_mean
        subtracted .+= c_mean
    end

    if time_mean
        α_approx = time_fit_iterate(X,GE.cell_type_sizes,GE.cell_type_labels,GE.embryo_id)
        X .-= α_approx
        subtracted .+= α_approx
    end

    if mother_mean
        ind_id = Int64.(floor.(GE.embryo_id/10))
        mother_mean = cluster_mean(X,ind_id) 
        X .-= mother_mean
        subtracted .+= mother_mean
    end
    return X, subtracted
end

function emb_shuffle_deep(GE::GeneExpressionStage;rem_cluster_mean=true,mother_mean=true,
        time_mean=true,shuff_vector=GE.cell_type_labels)
    # Shuffles embryos after time/mother/cell_type removal
    GE_copy = deepcopy(GE)
    X, subtracted = remove_variance(GE;rem_cluster_mean=rem_cluster_mean,mother_mean=mother_mean,time_mean=time_mean)

    e_bar = GE_copy.error_bars
    for i in unique(shuff_vector)
        f = findall(isequal(i),shuff_vector)
        sf = shuffle(f)
        X[:,f] .= X[:,sf]
        e_bar[:,f] .= e_bar[:,sf]
    end
    
    GE_copy.counts .= X .+ subtracted
    GE_copy.error_bars .= e_bar
    return GE_copy
end

function cluster_mean(GE::GeneExpressionStage)
    return hcat([mean(GE.counts[:,c .== GE.cell_type_labels],dims=2) for c in GE.cell_type_labels]...)
end

function cluster_mean(X,cell_type_labels)
    return hcat([mean(X[:,c .== cell_type_labels],dims=2) for c in cell_type_labels]...)
end


function time_fit_iterate(X,cell_type_sizes,cell_type_labels,embryo_id)
    nt = length(cell_type_sizes)
    α = randn(size(X,1),nt)
    α ./= sqrt(sum(α.^2))
    t = randn(length(unique(embryo_id)))

    e_id_dict = Dict(unique(embryo_id) .=> 1:length(unique(embryo_id)))
    idx_arr = [findall(isequal(i),cell_type_labels) for i in 1:nt]
    idx_emb = [findall(isequal(i),embryo_id) for i in unique(embryo_id)]
    α_vec = zeros(size(X))

    for _ = 1:100 # more iterations than necessary
        t_vec = [t[e_id_dict[e]] for e in embryo_id]
        for i = 1:nt
            α[:,i] .= (X[:,idx_arr[i]]*t_vec[idx_arr[i]])/sum(t_vec[idx_arr[i]].^2)
        end
        α ./= sqrt(sum(α.^2))

        α_vec .=  α[:,cell_type_labels]
        for i in 1:length(idx_emb)
            t[i] = sum(X[:,  idx_emb[i]].*α_vec[:,  idx_emb[i]]) / sum(α_vec[:,  idx_emb[i]].^2) 
        end
    end
    α_approx = [t[e_id_dict[e]] for e in embryo_id]' .*α[:,cell_type_labels]
    return α_approx
end
