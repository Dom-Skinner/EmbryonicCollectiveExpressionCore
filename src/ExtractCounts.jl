using CSV, DataFrames
using HDF5


function rough_norm(GE::GeneExpressionStage)
    count_raw = GE.counts    
    count_raw_indel = Float64.(copy(count_raw))
    count_raw_indel[count_raw_indel.==0] .= 0.5
    count_rough_norm = log.(count_raw_indel./sum(count_raw,dims=1) *1e6)# counts in RPM
    return GeneExpressionStage( GE.stage, "rough_norm", GE.genename, GE.embryo_id, count_rough_norm) 
end

function sanity_norm_core(df::DataFrame)
    bin_dir = (@__DIR__)*"/../bin/"
    
    g  = read(`uname`,String)
    if occursin("Linux",g)
        Sanity = bin_dir*"Sanity"
    else
        Sanity = bin_dir*"Sanity_macOS"
    end

    bin_rand = bin_dir*string(rand())[3:end]*"/" # to avoid race conditions (which happen otherwise when multiple snakemake jobs call this function)
    mkdir(bin_rand)
    raw_counts = bin_rand*"raw_counts.txt"
    output = bin_rand*"log_transcription_quotients.txt"
    err_bars = bin_rand*"ltq_error_bars.txt"
    
    CSV.write(raw_counts,df,delim="\t")
    run(`$Sanity -f $raw_counts -d $bin_rand`)
    df_norm = CSV.read(output,DataFrame,delim="\t")
    pred_error = CSV.read(err_bars,DataFrame,delim="\t")
    # clean up files
    rm(raw_counts,force=true)
    rm(output,force=true)
    rm(err_bars,force=true)
    rm(bin_rand)
    return df_norm, pred_error
end


function sanity_norm(GE::GeneExpressionStage)
    
    gene_dict = Dict("Gene".*string.(1:length(GE.genename)) .=> GE.genename)
    df = DataFrame("GeneID" => "Gene".*string.(1:length(GE.genename)))
    for i = 1:size(GE.counts,2)
        df[!,"Cell_"*string(i)] = Int64.(GE.counts[:,i])
    end
    
    df_norm,pred_error = sanity_norm_core(df)
    error_bars = Matrix(pred_error[:,2:end])
    gene_name = [gene_dict[x] for x in df_norm.GeneID]
    return GeneExpressionStage(GE.stage, "Sanity",gene_name, GE.embryo_id, Matrix(df_norm[:,2:end]),error_bars)
    
end

function sanity_norm_cluster(GE::GeneExpressionStage,clusters::Array)

    g = x-> string(split(x,"|")[2])
    df = DataFrame("GeneID" => g.(GE.genename))
    for i = 1:size(GE.counts,2)
        df[!,"Cell_"*string(i)] = Int64.(GE.counts[:,i])
    end

    df_total = []
    err_total = []
    for i in unique(clusters)
        df_r = df[:,vcat(1,findall(clusters.==i).+1)]
        df_norm,pred_error = sanity_norm_core(df_r)
        push!(df_total,df_norm)
        push!(err_total,pred_error)
    end

    keep_gene_names = intersect([d.GeneID for d in df_total]...)
    keep_genes = [findfirst(isequal(x), df.GeneID) for x in keep_gene_names]


    g_id = y ->[findfirst(isequal(x), y.GeneID) for x in keep_gene_names]
    counts = hcat([Matrix(s[g_id(s),2:end]) for s in df_total]...)
    error_bars = hcat([Matrix(s[g_id(s),2:end]) for s in err_total]...)

    embryo_id_sorted = GE.embryo_id[vcat([findall(clusters.==i) for i in unique(clusters)]...)]
    GE_tmp =  GeneExpressionStage(GE.stage, "Sanity_cluster", GE.genename[keep_genes], embryo_id_sorted, counts,error_bars)
    return GeneExpressionStage(GE_tmp,clusters[vcat([findall(clusters.==i) for i in unique(clusters)]...)])

end

function normalize_(GE::GeneExpressionStage,protocol::String,data_path::String,clusters::Array)
    @assert GE.normalization == "none"
    return sanity_norm_cluster(GE,clusters)
end

function normalize_(GE::GeneExpressionStage,protocol::String,data_path::String)
    @assert GE.normalization == "none"

    if protocol == "rough_norm"
        return rough_norm(GE)
    elseif protocol == "sanity"
        return sanity_norm(GE)
    end
    println("Warning: Unkown normalization protocol")
    return GE
end

function extract_gene_expression(data_path;stage=8,normalization = "none",clusters=missing)
    

    stages = h5read(data_path,"stages")
    genename = h5read(data_path,"genename")
    count_raw = h5read(data_path,"count_raw")
    embryo_id = h5read(data_path,"embryo_id")


    GE = GeneExpressionStage(stage,"none",genename,embryo_id[stages.==stage],count_raw[:,stages.==stage])

    if normalization != "none"
        if ismissing(clusters)
            return normalize_(GE,normalization,data_path)
        else
            @assert normalization == "Sanity_cluster"
            return normalize_(GE,normalization,data_path,clusters)
        end
    end
    return GE
end

function gene_restrict(GE::GeneExpressionStage,gene_list::Array;randomize=false)
    if !any(isnothing.(gene_list))
        if ismissing(GE.error_bars)
            return GeneExpressionStage(GE.stage,GE.normalization,GE.genename[gene_list],GE.embryo_id,GE.counts[gene_list,:])
        else
            if ismissing(GE.cell_type_labels)
                return GeneExpressionStage(GE.stage,GE.normalization,GE.genename[gene_list],GE.embryo_id,GE.counts[gene_list,:],GE.error_bars[gene_list,:])
            else
                return GeneExpressionStage(GE.stage,GE.normalization,GE.genename[gene_list],GE.embryo_id,
                            GE.counts[gene_list,:],GE.cell_type_labels, GE.cell_type_sizes,GE.error_bars[gene_list,:],
                            length(GE.genename[gene_list]),GE.N_samples)
            end
        end
    elseif randomize
        # if some genenames are missing, we may wish to replace with random data for clustering purposes.
        genename = [isnothing(g) ? "randomized" : GE.genename[g] for g in gene_list]
        counts = randn(length(gene_list),size(GE.counts,2))
        counts[ .!isnothing.(gene_list),:] .= GE.counts[filter(!isnothing,gene_list),:]
        
        err_bars = randn(length(gene_list),size(GE.counts,2))
        err_bars[ .!isnothing.(gene_list),:] .= GE.error_bars[filter(!isnothing,gene_list),:]
        return GeneExpressionStage(GE.stage,GE.normalization,genename,GE.embryo_id,counts,err_bars)
        
    end
    error()
end

function cell_restrict(GE::GeneExpressionStage,cell_list::Array;stage=0)
    
    s = (stage > 0) ? stage : GE.stage
    if ismissing(GE.error_bars)
        return GeneExpressionStage(s,GE.normalization,GE.genename,GE.embryo_id[cell_list],GE.counts[:,cell_list])
    else
        if ismissing(GE.cell_type_labels)
            return GeneExpressionStage(s,GE.normalization,GE.genename,GE.embryo_id[cell_list],GE.counts[:,cell_list],GE.error_bars[:,cell_list])
        else
            return GeneExpressionStage(s,GE.normalization,GE.genename,GE.embryo_id[cell_list],GE.counts[:,cell_list],GE.cell_type_labels[cell_list],GE.cell_type_sizes,GE.error_bars[:,cell_list],GE.ng,GE.N_samples)
        end
    end
end

function get_type(GE::GeneExpressionStage)
    @assert !ismissing(GE.cell_type_labels)
    if ismissing(GE.cell_type_sizes) 
        local_clusters = [GE.cell_type_labels[GE.embryo_id .== e] for e in unique(GE.embryo_id)]
        type_ =  sort(local_clusters[argmax(length.(local_clusters))])
        @assert length(type_) == GE.stage
        return type_
    else
        return vcat([repeat([i],GE.cell_type_sizes[i]) for i = 1:length(GE.cell_type_sizes)]...)
    end
end

function get_type(PC_GE::PCA_GeneExpression)
    return get_type(PC_GE.GE)
end
