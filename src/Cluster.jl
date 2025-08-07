using Graphs
using Random
using StatsBase:percentile

cfd = (x,M) -> sum(avech(M) .<= x) / (size(M,1)*(size(M,1)-1))*2
PAC = M -> cfd(0.9,M) - cfd(0.1,M)

function cluster_modal(GE,K,ntrials)
    GE_modal = deepcopy(GE)
    GE_modal.error_bars .*= 0.25 # Dial down uncertainty - find best fit clustering, we can assess the uncertainty later
    clst_vec = []
    size_vec = []
    M_score = []
    for i = 1:ntrials
        clusters,cluster_sizes = k_means_cluster(GE_modal,nclusters=K,randomize=true) 
        M,_ = cluster_to_concensus_matrix(clusters)
        push!(clst_vec,clusters)
        push!(size_vec,cluster_sizes)
        push!(M_score,M)
    end

    score = zeros(ntrials)
    for i = 1:ntrials 
        score[i] = sum([sum( M_score[i] .== M) for M in M_score])
    end
    return clst_vec[argmax(score)],size_vec[argmax(score)]
end

function iterate!(RC::RecursiveCluster;ntrials=4_000,regularization=0.02, 
        max_gene_num=5,max_clusters=3,PAC_tol=0.5,n_iters=200,n_shuff_trials=20,
        percentile_cutoff=10,modal_num=20)
    if !any(RC.is_active_node)
        return
    end
    idx = findfirst(RC.is_active_node)

    ## Do pSMD
    GE_local = cell_restrict(RC.GE_full,RC.partition_vector[idx],stage=RC.size_vector[idx])
    genes = pSMD(GE_local;ntrials=ntrials, regularization=regularization)
    top_genes = sortperm(genes)[end:-1:1]
    GE_local = gene_restrict(GE_local,top_genes[1:max_gene_num])
    println("Finished the SMD step")
    flush(stdout)

    PAC_score_and_GeneRes = [restrict_by_cluster_PAC(deepcopy(GE_local),K,n_iters=n_iters,subsample=0.9) for K = 2:minimum([max_clusters,(GE_local.stage>>1)])]
    PAC_score = [p[2] for p in PAC_score_and_GeneRes]
    K = argmin(PAC_score) + 1
    println("PAC scores:")
    println(PAC_score)
    # Now determine if we think cluster is real
    if minimum(PAC_score) < PAC_tol # Don't bother checking if PAC score is unambigouously small
        cluster_real = true
    else
        GE_shuff = deepcopy(PAC_score_and_GeneRes[K-1][1])
        shuffle_score = zeros(n_shuff_trials)
        for j = 1:length(shuffle_score)
            shuffle!(GE_shuff.embryo_id)
            shuffle_score[j] = concensus_M(GE_shuff,K,n_iters=n_iters,subsample=0.9) |> PAC
        end
        println("Shuffle scores")
        println(shuffle_score)
        if minimum(PAC_score) < percentile(shuffle_score,percentile_cutoff)
            println("Cluster is real")
            cluster_real = true
        else
            println("Cluster not real")
            cluster_real = false
        end
    end

    if cluster_real
      
        RC.is_active_node[idx] = false
        RC.GE_part[idx] = PAC_score_and_GeneRes[K-1][1]
        RC.PAC_score[idx] = minimum(PAC_score)

        #clusters,cluster_sizes = k_means_cluster(PAC_score_and_GeneRes[K-1][1],nclusters=K)  
        clusters,cluster_sizes = cluster_modal(PAC_score_and_GeneRes[K-1][1],K,modal_num)
      
        for i = 1:K
            add_vertex!(RC.g)
            add_edge!(RC.g,idx,nv(RC.g))
            push!(RC.is_active_node,cluster_sizes[i]>2)
            push!(RC.PAC_score,0.0)
            push!(RC.is_terminal,cluster_sizes[i]==2)
            push!(RC.partition_vector,RC.partition_vector[idx][clusters.==i])
            push!(RC.size_vector,cluster_sizes[i])
            push!(RC.GE_part,missing)
        end
    else
        RC.is_active_node[idx] = false
        RC.is_terminal[idx] = true
    end
    @everywhere (GC.gc();GC.gc())
	GC.gc()
    flush(stdout)
end

function merge_tree(RC,idx)
    downstream_nbhs = filter(x->x>idx,neighbors(RC.g,idx))
    downstream_nbhs = vcat(downstream_nbhs,[filter(x->x>idx,neighbors(RC.g,nn)) for nn in downstream_nbhs]...) |> unique    

    g_new,vmap = induced_subgraph(RC.g, setdiff(1:nv(RC.g),downstream_nbhs))
    idx_new = findfirst(isequal(idx),vmap)
    is_terminal = RC.is_terminal[vmap]
    is_terminal[idx_new]  = true
    is_active_node = RC.is_active_node[vmap]
    PAC_score = RC.PAC_score[vmap]
    PAC_score[idx_new]  = 0
    partition_vector = RC.partition_vector[vmap]
    size_vector = RC.size_vector[vmap]
    GE_part = RC.GE_part[vmap]
    GE_part[idx_new]  = missing
    return RecursiveCluster(g_new,is_terminal,is_active_node,PAC_score,partition_vector,
        size_vector,deepcopy(RC.GE_full),GE_part)
end

function iterate_with_guide!(RC::RecursiveCluster,RC_ref::RecursiveCluster)
    if !any(RC.is_active_node)
        return
    end

    
    idx_set = findall(RC.is_active_node) # potential nodes
    
    # remove any nodes that are already terminal
    for i in idx_set
        if RC_ref.is_terminal[i]
            RC.is_active_node[i] = false
            RC.is_terminal[i] = true
            return
        end
    end

    # choose next node to iterate based on reference
    idx = idx_set[argmin([minimum(filter(x->x>i,neighbors(RC_ref.g,i) )) for i in idx_set])]
    
    nb = filter(x->x>idx,neighbors(RC_ref.g,idx) )
    K = length(nb)
    
    if K > 0
        GE_local = cell_restrict(RC.GE_full,RC.partition_vector[idx],stage=RC.size_vector[idx])
        genes_keep = RC_ref.GE_part[idx].genename
        
        GE_local = gene_restrict(GE_local,[findfirst(isequal(g),GE_local.genename) for g in genes_keep],randomize=true)
        clusters,cluster_sizes = k_means_cluster(GE_local,nclusters=K,cluster_prior=RC_ref.size_vector[nb])


        ##########################################################################################
        # Swap clusters if same size and size > 2
        i = 1
        j = 2
        while i < length(cluster_sizes)
            while (j > i) && (j <= length(cluster_sizes))
                if (cluster_sizes[i] > 2) && (cluster_sizes[i] == cluster_sizes[j])
                    
                    c_mean = hcat([mean(GE_local.counts[:,c .== clusters],dims=2)[:,1] for c in [i,j]]...)
                    GE_tmp = gene_restrict(RC_ref.GE_full,[findfirst(isequal(g),RC_ref.GE_full.genename) for g in genes_keep],randomize=true)
                    G_mean = hcat([mean(GE_tmp.counts[:, RC_ref.partition_vector[n]],dims=2)[:,1] for n in nb[[i,j]]]...)        
                    opt1 = sum(( c_mean .- G_mean).^2)
                    opt2 = sum(( c_mean[:,[2,1]] .- G_mean).^2)
                    if opt2 < opt1
                        clusters_1 = clusters.== i
                        clusters_2 = clusters.== j
                        
                        clusters[clusters_1] .= j
                        clusters[clusters_2] .= i
                        i = 1
                        j = 1 # restart the loop
                    end
                end
                j += 1
            end
            i += 1
            j = i+1
        end
        ##########################################################################################
        @assert all(RC_ref.size_vector[nb] .== cluster_sizes)

        RC.is_active_node[idx] = false
        RC.GE_part[idx] = GE_local
        for i = 1:K
            add_vertex!(RC.g)
            add_edge!(RC.g,idx,nv(RC.g))
            push!(RC.is_active_node,cluster_sizes[i]>2)
            push!(RC.is_terminal,cluster_sizes[i]==2)
            push!(RC.partition_vector,RC.partition_vector[idx][clusters.==i])
            push!(RC.size_vector,cluster_sizes[i])
            push!(RC.GE_part,missing)
        end
    else
        RC.is_active_node[idx] = false
        RC.is_terminal[idx] = true
    end
end


function get_clusters(RC)
    @assert !any(RC.is_active_node)
    t_states = findall(RC.is_terminal)
    clusters = zeros(Int64,size(RC.GE_full.counts,2))
    for i in 1:length(t_states)
        clusters[RC.partition_vector[t_states[i]]] .= i
    end
    @assert all(clusters .> 0)
    cluster_sizes = RC.size_vector[t_states]
    return clusters,cluster_sizes
end

function minimal_cluster_info(RC)
    @assert !any(RC.is_active_node)
    t_states = findall(RC.is_terminal)
    clusters = zeros(Int64,size(RC.GE_full.counts,2))
    for i in 1:length(t_states)
        clusters[RC.partition_vector[t_states[i]]] .= i
    end
    @assert all(clusters .> 0)
    all_used_genes = vcat([s.genename for s in RC.GE_part[.!ismissing.(RC.GE_part)]]...) |> unique
    all_used_genes = filter(.!isequal("randomized"),all_used_genes)

    all_used_genes_idx = [findfirst(isequal(g),RC.GE_full.genename) for g in all_used_genes]
    GE_total = gene_restrict(RC.GE_full,all_used_genes_idx)
    cluster_sizes = RC.size_vector[t_states]
    GE_total = GeneExpressionStage(GE_total,clusters,cluster_sizes)
    return GE_total
end

function recluster(RC;randomize=true,subsample = 1.0)
    GE = minimal_cluster_info(RC)
    if randomize
        GE.counts .= GE.counts .+ GE.error_bars.*randn(size(GE.counts))
    end
    cells_to_keep = rand(size(GE.counts,2)) .< subsample
    GE = cell_restrict(GE,findall(cells_to_keep))

    RC_c = RecursiveCluster(GE)
    while any(RC_c.is_active_node)
        iterate_with_guide!(RC_c,RC)
    end
    
    clusters,sz = get_clusters(RC_c)
    clusters_all = zeros(Int64,length(cells_to_keep))
    clusters_all[cells_to_keep] .= clusters
    
    return clusters_all,sz
end

function concensus_M_block(RC::RecursiveCluster,subsample,block_len)
    n = size(RC.GE_full.counts,2)
    M_tot = zeros(n,n)
    I_tot = diagm(ones(n))
    for i = 1:block_len
            clusters,_ = recluster(RC;randomize=true,subsample = subsample)
            M,I = cluster_to_concensus_matrix(clusters)
            M_tot .+= M
            I_tot .+= I
    end
    return M_tot,I_tot
end


function concensus_M(RC::RecursiveCluster,nclusters;n_iters=1000,subsample=1.0)
    # for consistency with the GeneExpressionStage definition
    return concensus_M(RC,n_iters=n_iters,subsample=subsample)
end


function concensus_M(RC::RecursiveCluster;n_iters=1000,subsample=1.0)
    n = size(RC.GE_full.counts,2)
    block_len=1
    MI = pmap(x-> concensus_M_block(RC,subsample,block_len), 1:block_len:n_iters)
    @everywhere (GC.gc();GC.gc())
	GC.gc()
    flush(stdout)
    M_tot = sum([m[1] for m in MI])
    I_tot = sum([m[2] for m in MI]) + diagm(ones(n))
    return M_tot./I_tot
end
    


function part_way_cluster(RC,idx;randomize=true)
    # return the clustering for a partially completed recursive cluster
    c = zeros(Int64,length(RC.partition_vector[1]))
    nn = filter(x->x > idx,neighbors(RC.g,idx))
    for i = 1:length(nn)
        c[RC.partition_vector[nn[i]]] .= i
    end
    c = c[RC.partition_vector[idx]]
    if randomize
        X = RC.GE_part[idx].counts .+ RC.GE_part[idx].error_bars.*randn(size(RC.GE_part[idx].error_bars))
    else
        X = copy(RC.GE_part[idx].counts)
    end
    return X, c
end
