using Distributions
using Graphs
using MultivariateStats

struct CorrelationBounds
    C::Array
    C_cross::Array
    C_low::Array
    C_high::Array
    C_cross_low::Array
    C_cross_high::Array
end

struct MaxEntModel
    d::MvNormal
    nc::Int
    ng::Int
    genename::Array
    M::Array
    S::Array
    U::Array
    Ψ_sym::Array
    Ψu::Array
end

struct GeneExpressionStage
    stage::Int
    normalization::String
    genename::Array
    embryo_id::Array
    counts::Array
    cell_type_labels::Union{Array, Missing}
    cell_type_sizes::Union{Array, Missing}
    error_bars::Union{Array, Missing}
    ng::Int # Can be deduced from other elements, but included for convenience
    N_samples::Int
end



function GeneExpressionStage( stage::Int, normalization::String, genename::Array, embryo_id::Array, counts::Array)
    return GeneExpressionStage(stage, normalization, genename, embryo_id, counts,missing, missing,missing,size(counts,1),length(unique(embryo_id)))
end

function GeneExpressionStage( stage::Int, normalization::String, genename::Array, embryo_id::Array, counts::Array,error_bars::Array)
    return GeneExpressionStage(stage, normalization, genename, embryo_id, counts, missing,missing,error_bars,size(counts,1),length(unique(embryo_id)))
end

function GeneExpressionStage( GE::GeneExpressionStage, cell_type_labels::Array,cell_type_sizes::Array) # add cell type labels
    return GeneExpressionStage(GE.stage, GE.normalization, GE.genename, GE.embryo_id, GE.counts, cell_type_labels,cell_type_sizes,GE.error_bars,GE.ng,GE.N_samples)
end


struct RecursiveCluster
    g::SimpleGraph{T} where T
    is_terminal::BitVector
    is_active_node::BitVector
    PAC_score::Vector
    partition_vector::Vector{Vector{T}} where T
    size_vector::Vector
    GE_full::GeneExpressionStage
    GE_part::Vector{Union{GeneExpressionStage,Missing}}
end

function RecursiveCluster(GE_full)
    return RecursiveCluster(SimpleGraph(1),falses(1),trues(1),zeros(1),[Vector(1:size(GE_full.counts,2))],[GE_full.stage],GE_full,[missing])
end

function RecursiveCluster(GE_full,raw_data_path)
    # Restrict to only genes that have at least one read per embryo
    emb_id = h5read(raw_data_path,"embryo_id")
    stage_ID = h5read(raw_data_path,"stages")
    X = h5read(raw_data_path,"count_raw")
    genename = h5read(raw_data_path,"genename")
    gene_dict = Dict(genename .=> 1:length(genename))

    X = X[:,stage_ID.==GE_full.stage]
    emb_id = emb_id[stage_ID.==GE_full.stage]
    
    S = .*([sum(X[:,emb_id.==e],dims=2)[:,1] for e in unique(emb_id)]...) .>1
    GE_to_keep = gene_restrict(GE_full,findall(S[[gene_dict[g] for g in GE_full.genename]]))

    return RecursiveCluster(SimpleGraph(1),falses(1),trues(1),zeros(1),[Vector(1:size(GE_to_keep.counts,2))],[GE_to_keep.stage],GE_to_keep,[missing])
end

struct PCA_GeneExpression
    GE::GeneExpressionStage
    nPC::Int64
    M::PCA{T} where T
    rem_cluster_mean::Bool
    mother_mean::Bool
    time_mean::Bool
end

function PCA_GeneExpression(GE::GeneExpressionStage,nPC::Int64,rem_cluster_mean::Bool,mother_mean::Bool,time_mean::Bool)
    X,_ = remove_variance(GE,rem_cluster_mean=rem_cluster_mean,mother_mean=mother_mean,
        time_mean=time_mean,randomize=false)
    M = fit(PCA, X; maxoutdim=nPC);
    for i = 1:nPC # for consistency. For reasons that are beyond my faculties (numerical precision?) the PCA 
        # directions can be mirrored from one computer to another. Just fix them here, doesn't have to be a right handed coord system in any sense
        M.proj[:,i] .= M.proj[:,i]/sign(sum(M.proj[:,i])) 
    end
    return PCA_GeneExpression(GE,nPC,M,rem_cluster_mean,mother_mean,time_mean)
end
