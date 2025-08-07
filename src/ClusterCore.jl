using ScikitLearn
using Random
using StatsBase:mean
using HypothesisTests   
using Distributed     
using JuMP
import HiGHS
using Arpack
using Clustering
using Distances
#using MultivariateStats
using LinearAlgebra
using PyCall
const SVM = PyNULL()
function __init__()
    copy!(SVM, pyimport("sklearn.svm"))
end


# Need to run this to install code first time.
#    @sk_import svm: LinearSVC




function get_best_feature(model, X::Matrix, y,idx::BitVector)
    ScikitLearn.fit!(model, X[idx,:], y[idx])
 
    max_idx =  argmax(abs.(model.coef_)[1,:])
    return abs(model.coef_[1,max_idx]) > 1e-10 ? max_idx : 0
end


function run_trial!(g,regularization,X,clusters)
    idx = trues(size(X,1))
    y = copy(clusters)

    for i = 1:maximum(clusters), j = (i+1):maximum(clusters)
        idx .= (clusters .== i) .| (clusters .== j)
        model = SVM.LinearSVC(penalty="l1",dual=false,C=regularization,max_iter=2_000)
        feat = get_best_feature(model,X, y,idx)
        if feat > 0
            g[feat] += 1
        end
    end
end

function run_trial(regularization,X,clusters)
    g = zeros(size(X,2))
    run_trial!(g,regularization,X,clusters)
    return g
end

function cluster_assign(embryo_id,nc;cluster_merge = [])
    cluster_guess = zeros(Int64,length(embryo_id))
    cluster_assign!(cluster_guess,embryo_id,nc;cluster_merge = cluster_merge)
    return cluster_guess
end

function cluster_assign!(cluster_guess,embryo_id,nc;cluster_merge = [])
    
    if length(cluster_merge) == 0
        cluster_vector = repeat(1:(nc>>1),2)
    else
        @assert (sum(cluster_merge) << 1) == nc
        cluster_vector = vcat([repeat([i],cluster_merge[i]*2) for i = 1:length(cluster_merge)]...)
    end
    
    for em in unique(embryo_id)
        idx = findall(isequal(em),embryo_id)
        shuffle!(cluster_vector)
        cluster_guess[idx] .= cluster_vector[1:length(idx)]
    end
    
    return
end

function SMD(GE::GeneExpressionStage;ntrials=100, regularization=0.1,clust_fun! = ((c,GE) -> cluster_assign!(c,GE.embryo_id,GE.stage)) )
    
    X = permutedims(GE.counts)
    #emb_id = GE.embryo_id
    # First center the data. Should be approx variance stabilized already from normalization
    X_c = X .- mean(X,dims=1)
    
    g = zeros(size(X,2))
    
    clusters = ones(Int64,length(GE.embryo_id))

    #model = LinearSVC(solver=LIBSVM.Linearsolver.L1R_LR,cost=regularization) # instance
    
    for i = 1:ntrials
        # Take a new random cluster assignment
        clust_fun!(clusters,GE)
        # run trials to find best separating planes
        run_trial!(g,regularization,X_c,clusters)
    end
    
    return g
end

function pSMD(GE::GeneExpressionStage;ntrials=100, regularization=0.1,clust_fun = (GE -> cluster_assign(GE.embryo_id,GE.stage)) )
    
    g = zeros(size(GE.counts,1))
    
    mem_block = 200
    for i = 1:maximum([1,Int64(round(ntrials/mem_block))])
        clusters = [clust_fun(GE)  for i = 1:mem_block]
        X = permutedims(GE.counts)
        # First center the data. Should be approx variance stabilized already from normalization
        X_c = X .- mean(X,dims=1)
        g .+= sum(pmap(y->run_trial(regularization,X_c,y),clusters))
	@everywhere (GC.gc();GC.gc())
	GC.gc()
    flush(stdout)
    end
    
    return g
end

function transport_problem(PW,cluster_prior)
    # solve the linear transport problem to match 
    model = Model(HiGHS.Optimizer)        
    set_silent(model)
    @variable(model, x[1:size(PW,1),1:size(PW,2)] >= 0)        
    @objective(model,Min,sum(PW.*x))
    @constraint(model, [i in 1:size(PW,1)], sum(x[i, :]) == 1)
    @constraint(model, [j in 1:size(PW,2)], sum(x[:, j]) <= cluster_prior[j])
    optimize!(model)
    c = [argmax(value.(x[i,:])) for i = 1:size(PW,1)]
    return c
end

function iterated_k_means(X,nclusts)
    # run kmeans 100 times and pick the best one.
    score = [Inf]
    clst_best = zeros(Int64,size(X,2))
    for i = 1:100
        R = kmeans(X,nclusts; maxiter=500)
        clst = assignments(R)
        M = R.centers
        score_local = sum(x->x.^2,X .- M[:,clst])
        if score_local < score[1]
            score[1] = score_local
            clst_best .= clst
        end
    end
    return clst_best
end

function cluster_initial(X,emb_id, cluster_prior)
    nclusts = length(cluster_prior)
    
    #R = kmeans(X, nclusts; maxiter=500)
    #clst = assignments(R)
    clst = iterated_k_means(X,nclusts)
    
    #= hclust
    R = pairwise(Euclidean(),X , dims=2)
    hclu = hclust(R,linkage=:ward)
    clst = cutree(hclu,k=nclusts)
    =#
    # sort so that the first cluster has the fewest points, 2nd has the 2nd fewest etc.
    to_sort = sortperm([sum(clst .==k) for k in 1:nclusts])
    clst = invperm(to_sort)[clst]

    c_cents = hcat([mean(X[:,clst .==k],dims=2) for k = 1:nclusts]...)
    
    clusters= zeros(Int64,size(X,2))
    
    for id in unique(emb_id)
        PW = pairwise(Euclidean(),X[:,emb_id .==id],c_cents)
        clusters[emb_id .==id] .= transport_problem(PW,cluster_prior)
    end

    return clusters
end


function merge_cluster(X,clusters,emb_id,clst_size;cluster_prior=[])
    
    nclusters= length(clst_size)
    @assert nclusters == length(unique(clusters))
    d = distance_ward(X,copy(clusters),nclusters)
    
    while minimum(d) < Inf
        idx = argmin(d)
        cluster_size = copy(clst_size)
        c = copy(clusters)
      
        # merge
        cluster_size[idx[1]] += cluster_size[idx[2]]
        cluster_size = cluster_size[setdiff(1:nclusters,idx[2])]
        to_sort = sortperm(cluster_size)
        cluster_size = cluster_size[to_sort]
        
        c[c .== idx[2]] .= idx[1]
        c[c .> idx[2]] .-= 1
        c .= invperm(to_sort)[c]
        if (length(cluster_prior) == 0) || dynamic_feasibility(cluster_size,cluster_prior)
            c_final = k_means_HW_core(X,emb_id,cluster_size,c)
            return c_final, cluster_size
        else
            d[idx] = Inf
        end
    end
end

function distance_ward(X,clusters,nc)
    
    clusts = [X[:,c.==clusters] for c in 1:nc]
    f = x -> sum((x .- mean(x,dims=2)).^2)
    d = f.(clusts)
    d_ward = fill(Inf,nc,nc)
    for i = 1:nc
        for j = (i+1):nc
            d_ward[i,j] = f(hcat(clusts[i],clusts[j])) - d[i] - d[j]
        end
    end
    return d_ward
end

function φ(Sn)
    μ = mean(Sn,dims=2)
    return sum( (Sn .- μ).^2)
end


function construct_similarity_matrix(X,emb_id)
    D = pairwise(Euclidean(),X,dims=2)
    A = zeros(size(D))

    for idx = 1:size(A,1)
        for e in unique(emb_id)
            cond = findall(emb_id .== e)
            sp = sortperm(D[cond,idx])
            if length(sp) == 1
                continue
            end
            σ = (D[cond[sp[1]],idx] == 0) ? D[cond[sp[2]],idx] : D[cond[sp[1]],idx]
            for i = 1:(length(cond)>>1)
                A[cond[sp[(2*i-1):2*i]],idx] .= exp.(- 0.5*i*D[cond[sp[(2*i-1):2*i]],idx].^2/ σ^2)
            end
        end
        A[idx,idx] = 0
    end

    for i = 1:size(A,1)
        for j = (i+1):size(A,2)
            A[i,j] = maximum([A[i,j],A[j,i]])
            A[j,i] = A[i,j]
        end
    end

    Deg = [sum(A[i,:]) for i in 1:size(A,1)]
    A = diagm(1 ./sqrt.(Deg)) * A * diagm(1 ./sqrt.(Deg))
    return Symmetric(A)
end    

using Arpack

function merge_vec(vec,i,j)
    @assert i < j
    new_vec = vec[setdiff(1:length(vec),i)]
    new_vec[j-1] += vec[i]
    sort!(new_vec)
    return new_vec
end

function dynamic_feasibility(vec1,vec2)
    d = Dict()
    dynamic_feasibility!(d,vec1,vec2)
    return d[vec1]
end
function dynamic_feasibility!(d,vec1,vec2)
    if haskey(d,vec1)
        return
    end
    if length(vec1) == length(vec2)
        d[vec1] = all(vec1 .== vec2)
        return
    end

    for i = 1:length(vec1)
        for j = (i+1):length(vec1)
            vec3 = merge_vec(vec1,i,j)
            dynamic_feasibility!(d,vec3,vec2)
            if d[vec3]
                d[vec1] = true
                return
            end
        end
    end

    d[vec1] = false
    return
end

function k_means_cluster(GE;randomize=false,subsample=1.0,nclusters=0,stage=0,cluster_prior=[])
    
    
    stage = (stage > 0) ? stage : GE.stage
    nclusters = (nclusters > 0) ? nclusters : (stage >> 1)
    
    X  = copy(GE.counts)
    if randomize
        X .+= (GE.error_bars.*randn(size(GE.counts)))
    end
    cells_to_keep = rand(size(GE.counts,2)) .< subsample
    emb_id = GE.embryo_id[cells_to_keep]
    X = X[:,cells_to_keep]

    X_spc = copy(X)

    cluster_merge = 2*ones(Int64,(stage >> 1))

    clusters_init = cluster_initial(X_spc,emb_id, cluster_merge)
    
    clusters =  k_means_HW_core(X_spc,emb_id, cluster_merge,clusters_init)
    
    for i = 1:((stage>>1) - nclusters)
        clusters, cluster_merge = merge_cluster(X_spc,clusters,emb_id,cluster_merge,cluster_prior=cluster_prior)
    end

    clusters_all = zeros(Int64,size(GE.counts,2))
    clusters_all[cells_to_keep] .= clusters

    return clusters_all, cluster_merge
end

  
function k_means_HW(GE::GeneExpressionStage;randomize=false,subsample=1.0,nclusters=0)
    
    nclusters = (nclusters > 0) ? nclusters : (GE.stage >> 1)
    X  = copy(GE.counts)
    if randomize
        X .+= (GE.error_bars.*randn(size(GE.counts)))
    end

    cells_to_keep = rand(size(GE.counts,2)) .< subsample
    emb_id = GE.embryo_id[cells_to_keep]
    X = X[:,cells_to_keep]
    # K-means with the Hartigan-Wong method, except also ensuring that 
    # each cluster contains cells from each embryo.

    
    cluster_merge = 2*ones(Int64,(GE.stage >> 1))

    clusters_init = cluster_initial(X,emb_id, cluster_merge)
    
    clusters =  k_means_HW_core(X,emb_id, cluster_merge,clusters_init)
    
    for i = 1:((GE.stage>>1) - nclusters)
        clusters, cluster_merge = merge_cluster(X,clusters,emb_id,cluster_merge)
    end

    clusters_all = zeros(Int64,size(GE.counts,2))
    clusters_all[cells_to_keep] .= clusters

    return clusters_all
end

function k_means_HW_core(X,emb_id,cluster_true_size,clusters_init)
    
    clusters = copy(clusters_init)
    clust_num = length(cluster_true_size)
    emb_lookup = [findall(emb_id .== i) for i in unique(emb_id)]
    
    iters = [(i,n,m) for i in 1:length(emb_lookup), n in 1:clust_num for m in (n+1):clust_num]
    
    block_len = 1500
    for j = 1:500 # arbitrary number of iterations
        cmat = [findall(clusters.==n) for n = 1:clust_num]
        Δ_max = 0
        to_switch = [0,0,0,0]
        # loop through all possible switches and pick the one that improves the objective the most
        shuffle!(iters)
        
        for kk in iters[1:minimum([block_len,length(iters)])]#i in 1:length(emb_lookup), n in 1:clust_num, m in (n+1):clust_num
            i,n,m = kk
            cn = cmat[n]
            cm = cmat[m]#findall(clusters.==m)
            x_options = intersect(emb_lookup[i],cn)
            y_options = intersect(emb_lookup[i],cm)

            # First do swaps for x and y
            for x in x_options, y in y_options 
                Δ = φ(X[:,cn]) + φ(X[:,cm]) - φ(X[:,vcat(setdiff(cn,x),y)])- φ(X[:,vcat(setdiff(cm,y),x)])

                if Δ > Δ_max
                    Δ_max = Δ
                    to_switch = [n,m,x,y]
                end
            end
            
            # Now do moves from one cluster to another if clusters are unevenly sized

            if (length(y_options) < cluster_true_size[m]) 
                for z in x_options
                    Δ = φ(X[:,cn]) + φ(X[:,cm]) - φ(X[:,vcat(cm,z)])- φ(X[:,setdiff(cn,z)])
                    if Δ > Δ_max
                        Δ_max = Δ
                        to_switch = [m,m,z,z] # z was in n now in m
                    end
                end
            end


            if (length(x_options) < cluster_true_size[n]) 
                for z in y_options
                    Δ = φ(X[:,cn]) + φ(X[:,cm]) - φ(X[:,vcat(cn,z)])- φ(X[:,setdiff(cm,z)])
                    if Δ > Δ_max
                        Δ_max = Δ
                        to_switch = [n,n,z,z] # z was in n now in m
                    end
                end
            end

        end


        if Δ_max == 0
            if block_len >= length(iters)
                return clusters
            else
                block_len = block_len*2
            end
        else
            # make the switch
            n,m,x,y = to_switch
            clusters[y] = n
            clusters[x] = m
        end
        
    end
    println("Warning: did not converge")
    return clusters
end

function concensus_M_block(GE::GeneExpressionStage,subsample,nclusters,block_len;stage=stage)
    n = size(GE.counts,2)
    M_tot = zeros(n,n)
    I_tot = diagm(ones(n))
    for i = 1:block_len
            clusters,_ = k_means_cluster(GE,randomize=true,subsample=subsample,nclusters=nclusters,stage=stage)
            M,I = cluster_to_concensus_matrix(clusters)
            M_tot .+= M
            I_tot .+= I
    end
    return M_tot,I_tot
end
function concensus_M(GE::GeneExpressionStage,nclusters;n_iters=1000,subsample=1.0,stage = 0)
    n = size(GE.counts,2)
    block_len=1
    MI = pmap(x-> concensus_M_block(GE,subsample,nclusters,block_len,stage=stage), 1:block_len:n_iters)
    @everywhere (GC.gc();GC.gc())
	GC.gc()
    flush(stdout)
    M_tot = sum([m[1] for m in MI])
    I_tot = sum([m[2] for m in MI]) + diagm(ones(n))
  
    return M_tot./I_tot
end
    

function cluster_to_concensus_matrix(clusters)
    n = length(clusters)
    M = zeros(n,n)
    I = zeros(n,n)
    for i = 1:n, j = (i+1):n
            if (clusters[i] != 0) .& (clusters[j] != 0)
                    I[i,j] = 1
                    I[j,i] = 1
            end
            if (clusters[i] != 0) .& (clusters[i] == clusters[j])                        
                    M[i,j] = 1
                    M[j,i] = 1
            end
    end
    return M,I
end

function restrict_by_cluster(GE,nclusters;n_iters=50,subsample=1.0,stage=0)
    # Assign an overall cluster consistency score to a set of genes, 
    # and see if removing any genes decreases the score (more consistent).
    if length(GE.genename) == 1
        return GE
    end
    stage = (stage > 0) ? stage : GE.stage
    M_score = x -> sum( x .*(1 .-x))        
    ng = length(GE.genename)
    println("There are ", ng, " genes remaining")
    M0 = concensus_M(GE,nclusters,n_iters=n_iters,subsample=subsample,stage=stage) |> M_score
    M_mat = [ concensus_M(gene_restrict(GE,setdiff(1:ng,i)),nclusters,n_iters=n_iters,subsample=subsample,stage=stage) |> M_score for i = 1:ng]
    if minimum(M_mat) <= M0*(1 + sqrt(n_iters)/n_iters) # allow some tolerance but seek to get sparser
            return restrict_by_cluster(gene_restrict(GE,setdiff(1:ng,argmin(M_mat))),nclusters,n_iters=n_iters,subsample=subsample,stage=stage)
    else
            return GE
    end
end



function restrict_by_cluster_PAC(GE,nclusters;n_iters=50,subsample=1.0,stage=0)
    # Assign an overall cluster consistency score to a set of genes, 
    # and see if removing any genes decreases the score (more consistent).
   
    stage = (stage > 0) ? stage : GE.stage

    cfd = (x,M) -> sum(avech(M) .<= x) / (size(M,1)*(size(M,1)-1))*2
    PAC = M -> cfd(0.9,M) - cfd(0.1,M)
    ng = length(GE.genename)
    println("There are ", ng, " genes remaining")
    
    if length(GE.genename) == 1
        M0 = concensus_M(GE,nclusters,n_iters=n_iters,subsample=1.0,stage=stage) |> PAC
        return GE, M0
    end
    M0 = concensus_M(GE,nclusters,n_iters=n_iters,subsample=subsample,stage=stage) |> PAC
    M_mat = [ concensus_M(gene_restrict(GE,setdiff(1:ng,i)),nclusters,n_iters=n_iters,subsample=subsample,stage=stage) |> PAC for i = 1:ng]
    if minimum(M_mat) <= M0*(1 + 2*sqrt(n_iters)/n_iters) # allow some tolerance but seek to get sparser
        println("One gene removed")
            return restrict_by_cluster_PAC(gene_restrict(GE,setdiff(1:ng,argmin(M_mat))),nclusters,n_iters=n_iters,subsample=subsample,stage=stage)
    else
        if length(GE.genename) > 2
            gene_importance = sortperm(M_mat)
            println("multiple genes removed")
            M_mat_red = [ concensus_M(gene_restrict(GE,setdiff(1:ng,gene_importance[1:i])),nclusters,n_iters=n_iters,subsample=subsample,stage=stage) |> PAC for i = 2:(length(gene_importance)-1)]
            if minimum(M_mat_red) <= M0*(1 + 2*sqrt(n_iters)/n_iters) # allow some tolerance but seek to get sparser
                return restrict_by_cluster_PAC(gene_restrict(GE,setdiff(1:ng,gene_importance[1:(argmin(M_mat_red)+1)])),nclusters,n_iters=n_iters,subsample=subsample,stage=stage)
            end
        end

        M0 = concensus_M(GE,nclusters,n_iters=n_iters,subsample=1.0,stage=stage) |> PAC
        return GE, M0   
    end
end

function restrict_by_cluster_iter(GE,nclusters;n_iters=50,subsample=1.0)
    # Assign an overall cluster consistency score to a set of genes, 
    # and see if removing any genes decreases the score (more consistent).
    M_score = x -> sum( x .*(1 .-x))        
    ng = length(GE.genename)
    println("There are ", ng, " genes remaining")
    M0 = concensus_M(GE,nclusters,n_iters=n_iters) |> M_score
    for i = ng:-1:1
        M_score_new = concensus_M(gene_restrict(GE,setdiff(1:ng,i)),nclusters,n_iters=n_iters,subsample=subsample) |> M_score
    
        if M_score_new <= M0
                return restrict_by_cluster_iter(gene_restrict(GE,setdiff(1:ng,i)),nclusters,n_iters=n_iters,subsample=subsample)
        end
    end
    return GE
end

function cluster_DE_test(GE)
    @assert !ismissing(GE.cell_type_labels)
    return cluster_DE_test(GE,GE.cell_type_labels) 
end
function cluster_DE_test(GE,clusters)
    # Find genes that are differentially expressed across clusters
    cnum = length(unique(clusters))
    cnum_pair =  (cnum*(cnum+1)) >> 1
    p_val_arr = zeros(length(GE.genename))
    for gene = 1:length(GE.genename)
            p_arr = []
            for c1 in clusters
                    for c2 in clusters
                        if c1 < c2
                            X1 = GE.counts[gene,clusters.==c1]        
                            X2 = GE.counts[gene,clusters.==c2]     
                            p = MannWhitneyUTest(X1,X2) |> pvalue 
                            push!(p_arr,cnum_pair* p*length(GE.genename))
                        end
                    end
            end
            p_val_arr[gene] = minimum(p_arr)
    end
    return p_val_arr
end

function maternal_test(GE,clusters,germ_cell)
    # Find genes that are differentially expressed across clusters
    p_val_arr = zeros(length(GE.genename))
    X = copy(GE.counts)
    embryo_mean = cluster_mean(X,GE.embryo_id) 
    X .-= embryo_mean

    for gene = 1:length(GE.genename)
        X1 = X[gene,clusters.==germ_cell]        
        X2 = X[gene,clusters.!=germ_cell]     
        p_val_arr[gene] = MannWhitneyUTest(X1,X2) |> pvalue 
    end
    return p_val_arr
end
