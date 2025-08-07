using HDF5

function GE_write_convenience(file,GE,str="")
    write(file, str*"Type", "GeneExpressionStage")
    write(file, str*"stage", GE.stage)
    write(file, str*"normalization", GE.normalization)
    write(file, str*"embryo_id", GE.embryo_id)
    write(file, str*"genename", String.(GE.genename))
    write(file, str*"counts", GE.counts)
    write(file, str*"ng", GE.ng)
    write(file, str*"N_samples", GE.N_samples)
    if ismissing(GE.cell_type_labels)
        write(file, str*"cell_type_labels", "missing")
    else
        write(file, str*"cell_type_labels", GE.cell_type_labels)
    end
    if ismissing(GE.cell_type_sizes)
        write(file, str*"cell_type_sizes", "missing")
    else
        write(file, str*"cell_type_sizes", GE.cell_type_sizes)
    end

    if ismissing(GE.error_bars)
        write(file, str*"error_bars", "missing")
    else
        write(file, str*"error_bars", GE.error_bars)
    end
end

function save(save_str,GE::GeneExpressionStage)
    h5open(save_str, "w") do file
        GE_write_convenience(file,GE)
    end
end

function save(save_str,Mod::MaxEntModel)
    h5open(save_str, "w") do file
        write(file, "Type", "MaxEntModel")
        write(file, "stage", Mod.nc)
        write(file, "ng", Mod.ng)
        write(file, "mean", Mod.d.μ)
        write(file, "covariance", Matrix(Mod.d.Σ))
        write(file, "genename", String.(Mod.genename))
        write(file, "M", vcat(Mod.M...))
        
        if length(Mod.Ψ_sym) > 0
            write(file, "S", vcat(Mod.S...))
            write(file, "Psi_sym", vcat(Mod.Ψ_sym...))
            write(file, "Psi_sym_decode", size.(Mod.Ψ_sym,1))
        end
        if length(Mod.Ψu) > 0
            write(file, "U", vcat(Mod.U...))
            write(file, "Psi_U", vcat(Mod.Ψu...))
            write(file, "Psi_U_decode", size.(Mod.Ψu,1))
        end
        
    end
end

function save(save_str,C_tot::CorrelationBounds)
    h5open(save_str, "w") do file
        write(file, "Type", "CorrelationBounds")
        write(file, "C", vcat(C_tot.C...))
        write(file, "C_cross", vcat(C_tot.C_cross...))
        write(file, "C_low", vcat(C_tot.C_low...))
        write(file, "C_cross_low", vcat(C_tot.C_cross_low...))
        write(file, "C_high", vcat(C_tot.C_high...))
        write(file, "C_cross_high", vcat(C_tot.C_cross_high...))
    end
end

function save(save_str,RC::RecursiveCluster)
    h5open(save_str, "w") do file
        write(file, "Type", "RecursiveCluster")
        write(file, "Graph",hcat(src.(collect(edges(RC.g))),dst.(collect(edges(RC.g)))))
        write(file, "is_terminal",Vector(RC.is_terminal))
        write(file, "is_active_node",Vector(RC.is_active_node))
        write(file, "PAC_score",RC.PAC_score)
        write(file, "partition_vector",vcat(RC.partition_vector...))
        write(file, "partition_vector_length",length.(RC.partition_vector))
        write(file, "size_vector",RC.size_vector)
        GE_write_convenience(file,RC.GE_full,"GE_full_")
        for  i = 1:length(RC.GE_part)
            if ismissing(RC.GE_part[i])
                write(file,"GE_part_"*string(i)*"_Type","missing")
            else
                GE_write_convenience(file,RC.GE_part[i],"GE_part_"*string(i)*"_")
            end
        end
    end
end

function load(save_str)
    input_type = h5read(save_str,"Type")
    if input_type == "GeneExpressionStage"
        return load_GeneExpressionStage(save_str)
    elseif input_type == "MaxEntModel"
        return load_MaxEntModel(save_str)
    elseif input_type == "RecursiveCluster"
        return load_RecursiveCluster(save_str)
    elseif input_type == "CorrelationBounds"
        return load_CorrelationBounds(save_str)
    end
end

function GE_read_convenience(save_str,str="")
    stage = h5read(save_str,str*"stage")
    normalization = h5read(save_str,str*"normalization")
    embryo_id = h5read(save_str,str*"embryo_id")
    counts = h5read(save_str,str*"counts")
    genename = h5read(save_str,str*"genename")
    ng = h5read(save_str,str*"ng")
    N_samples = h5read(save_str,str*"N_samples")
    cell_type_labels = h5read(save_str,str*"cell_type_labels")
    try
        cell_type_sizes = h5read(save_str,str*"cell_type_sizes")
    catch
        cell_type_sizes = "missing"
    end
    error_bars = h5read(save_str,str*"error_bars")
    
    if cell_type_labels == "missing"
        cell_type_labels = missing
    end

    if cell_type_sizes == "missing"
        cell_type_sizes = missing
    end

    if error_bars == "missing"
        error_bars = missing
    end

    return GeneExpressionStage(stage, normalization, genename, embryo_id, counts, 
        cell_type_labels,cell_type_sizes, error_bars,ng,N_samples)
end
function load_GeneExpressionStage(save_str)
    return GE_read_convenience(save_str)
    stage = h5read(save_str,"stage")
    normalization = h5read(save_str,"normalization")
    embryo_id = h5read(save_str,"embryo_id")
    counts = h5read(save_str,"counts")
    genename = h5read(save_str,"genename")
    ng = h5read(save_str,"ng")
    N_samples = h5read(save_str,"N_samples")
    cell_type_labels = h5read(save_str,"cell_type_labels")
    try
        cell_type_sizes = h5read(save_str,"cell_type_sizes")
    catch
        cell_type_sizes = "missing"
    end
    error_bars = h5read(save_str,"error_bars")
    
    if cell_type_labels == "missing"
        cell_type_labels = missing
    end

    if cell_type_sizes == "missing"
        cell_type_sizes = missing
    end

    if error_bars == "missing"
        error_bars = missing
    end

    return GeneExpressionStage(stage, normalization, genename, embryo_id, counts, 
        cell_type_labels,cell_type_sizes, error_bars,ng,N_samples)
end



function load_CorrelationBounds(save_str)
    quick_fun = (C,ng) ->[C[i:i+(ng-1),:] for i = 1:ng:(size(C,1)-ng+1)]
    C = h5read(save_str,"C")
    ng = size(C,2)
    C = quick_fun(C,ng)
    C_low = quick_fun(h5read(save_str,"C_low"),ng)
    C_high = quick_fun(h5read(save_str,"C_high"),ng)
    C_cross = quick_fun(h5read(save_str,"C_cross"),ng)
    C_cross_low = quick_fun(h5read(save_str,"C_cross_low"),ng)
    C_cross_high = quick_fun(h5read(save_str,"C_cross_high"),ng)
    return CorrelationBounds(C,C_cross,C_low,C_high,C_cross_low,C_cross_high)
end

function load_MaxEntModel(save_str)

    
    nc = h5read(save_str, "stage")
    ng = h5read(save_str, "ng")
    μ = h5read(save_str, "mean")
    Σ = h5read(save_str, "covariance")
    genename = h5read(save_str, "genename")
    
    M = h5read(save_str, "M")
    nt = Int(size(M,1)/size(M,2))
    M = [M[ (1+(i-1)*ng): i*ng,:] for i = 1:nt]

    
    has_key = true
    try
        S = h5read(save_str, "S")
    catch
       has_key = false
    end
    S = has_key ? h5read(save_str, "S") : []
    Ψ_sym_decode = has_key ? h5read(save_str, "Psi_sym_decode") : []
    Ψ_sym = has_key ? h5read(save_str, "Psi_sym") : []
    
    
    S = [S[ (1+(i-1)*ng): i*ng,:] for i = 1:length(Ψ_sym_decode)]
    Ψ_sym = [Ψ_sym[ (i == 0 ? 1 : (sum(Ψ_sym_decode[1:i])+1)):sum(Ψ_sym_decode[1:i+1]),:] for i = 0:(length(Ψ_sym_decode)-1)]

    has_key = true
    try
        U = h5read(save_str, "U")
    catch
        has_key = false
    end
    U = has_key ? h5read(save_str, "U") : []
    Ψ_U = has_key ? h5read(save_str, "Psi_U") : []
    Ψ_U_decode = has_key ? h5read(save_str, "Psi_U_decode") : []

    U = [U[ (1+(i-1)*ng): i*ng,:] for i = 1:length(Ψ_U_decode)]
    Ψ_U = [Ψ_U[ (i == 0 ? 1 : (sum(Ψ_U_decode[1:i])+1)):sum(Ψ_U_decode[1:i+1]),:] for i = 0:(length(Ψ_U_decode)-1)]

    return MaxEntModel(MvNormal(μ, Σ), nc,ng,genename,M,S,U,Ψ_sym,Ψ_U)
end


function load_RecursiveCluster(save_str)
    E = h5read(save_str, "Graph")
    g = graph_construct(E)
    is_terminal = BitVector(h5read(save_str, "is_terminal"))
    is_active_node = BitVector(h5read(save_str, "is_active_node"))
    PAC_score = zeros(length(is_terminal))
    try
        PAC_score .= h5read(save_str, "PAC_score")    
    catch
        println("no PAC defined yet")
    end
    

    partition_vector = h5read(save_str, "partition_vector")
    partition_vector_length = h5read(save_str, "partition_vector_length")
    pv = vcat([0],cumsum(partition_vector_length))
    partition_vector = [partition_vector[pv[i-1]+1:pv[i]] for i = 2:length(pv)]

    size_vector = h5read(save_str, "size_vector")
    
    GE_full = GE_read_convenience(save_str,"GE_full_")
    GE_part = Vector{Union{GeneExpressionStage,Missing}}(undef,length(size_vector))
        
    
    for  i = 1:length(GE_part)
        m  = h5read(save_str,"GE_part_"*string(i)*"_Type")
        if m == "missing"
            GE_part[i] = missing
        else
            GE_part[i] = GE_read_convenience(save_str,"GE_part_"*string(i)*"_")
        end
    end
    
    return RecursiveCluster(g,is_terminal,is_active_node,PAC_score,partition_vector,size_vector,GE_full,GE_part)
end