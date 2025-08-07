module EmbryonicCollectiveExpressionCore

include("ClassDefinitions.jl")
include("ExtractCounts.jl")
include("Utils.jl")
include("EmpiricalStats.jl")
include("ClusterCore.jl")
include("Cluster.jl")
include("MaxEnt.jl")
include("Ranking.jl")
include("ContactGraph.jl")
include("ReadWrite.jl")

export
    # ExtractCounts
    GeneExpressionStage,extract_gene_expression, 
    gene_restrict, cell_restrict, get_type, sanity_norm_core, PCA_GeneExpression,
    # Utils
    vech, avech, ivech, iavech, pos_def_cone, indicator,
    triangle_inv, arrange_GE_data, setdiff_with_repeats,generate_synthetic_data,embryo_bootstrap,
    emb_shuffle_deep,cluster_mean, permute_mat,
    # EmpiricalStats
    bootstrap_correlation, CorrelationBounds, empirical_mean_corr, exact_4_point, 
    empirical_4_point, bootstrap_4_correlation, Gaussian_4_point,
    # SMD
    SMD,pSMD,k_means_HW,k_means_HW_core,cluster_initial,restrict_by_cluster,concensus_M_block, recluster,
    restrict_by_cluster_iter,cluster_DE_test,maternal_test,concensus_M,merge_cluster,k_means_cluster,
    restrict_by_cluster_PAC,
    # Cluster
    RecursiveCluster,iterate!,minimal_cluster_info,iterate_with_guide!,get_clusters,
    cluster_to_concensus_matrix,merge_tree,part_way_cluster,
    # MaxEnt
    MaxEntModel, rand_sample, fit_independent_max_ent,
    exact_mean_corr, a_adj_mat, adj_mat, graph_construct, pos_def_proj_Γij, 
    params_to_model,obj_fun_SU, pos_def_proj_Γij!,exact_within_corr, max_ent_restricted,
    find_maxent_params_SU,
    # Ranking
    CV_array_to_rank,
    # ContactGraph
    index_to_cell_name, get_contact_graph,print_gene_ranking,
    # ReadWrite
    save,load
end 
