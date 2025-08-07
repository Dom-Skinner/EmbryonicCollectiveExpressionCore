# This file contains hard coded information about the contact graph

function index_to_cell_name_64(i)
    arr = Array{String}(undef,64)
    arr[1] =  arr[2] = "B7.6"
    arr[3] =  arr[4] = "A7.8"
    arr[5] =  arr[6] = "B7.5"
    arr[7] =  arr[8] = "B7.2"
    arr[9] =  arr[10]= "B7.1"
    arr[11] = arr[12]= "B7.8"
    arr[13] = arr[14]= "b7.14"
    arr[15] = arr[16]= "b7.13"
    arr[17] = arr[18]= "b7.16"
    arr[19] = arr[20]= "b7.11"
    arr[21] = arr[22]= "b7.15"
    arr[23] = arr[24]= "b7.12"
    arr[25] = arr[26]= "B7.7"
    arr[27] = arr[28]= "B7.4"
    arr[29] = arr[30]= "B7.3"
    arr[31] = arr[32]= "b7.9"
    arr[33] = arr[34]= "b7.10"
    arr[35] = arr[36]= "a7.16"
    arr[37] = arr[38]= "a7.15"
    arr[39] = arr[40]= "a7.14"
    arr[41] = arr[42]= "a7.12"
    arr[43] = arr[44]= "a7.11"
    arr[45] = arr[46]= "a7.13"
    arr[47] = arr[48]= "A7.6"
    arr[49] = arr[50]= "A7.1"
    arr[51] = arr[52]= "A7.5"
    arr[53] = arr[54]= "A7.2"
    arr[55] = arr[56]= "A7.7"
    arr[57] = arr[58]= "A7.3"
    arr[59] = arr[60] = "a7.10"
    arr[61] = arr[62] = "a7.9"
    arr[63] = arr[64] = "A7.4"
    return arr[i]
end


function index_to_cell_name_32(i)
    arr = Array{String}(undef,32)
    arr[1] = arr[2] = "A6.1"
    arr[3] = arr[4] = "A6.3"
    arr[5] = arr[6] = "A6.2"
    arr[7] = arr[8] = "A6.4"
    arr[9] = arr[10] = "a6.5"
    arr[11] = arr[12] = "a6.6"
    arr[13] = arr[14] = "a6.7"
    arr[15] = arr[16] = "a6.8"
    arr[17] = arr[18] = "b6.5"
    arr[19] = arr[20] = "b6.6"
    arr[21] = arr[22] = "b6.8"
    arr[23] = arr[24] = "b6.7"
    arr[25] = arr[26] = "B6.2"
    arr[27] = arr[28] = "B6.4"
    arr[29] = arr[30] = "B6.1"
    arr[31] = arr[32] = "B6.3"
    return arr[i]
end


function index_to_cell_name_16(i)
    arr = Array{String}(undef,16)
    arr[1] = arr[2] = "A5.1"
    arr[3] = arr[4] = "A5.2"
    arr[5] = arr[6] = "a5.3"
    arr[7] = arr[8] = "a5.4"
    arr[9] = arr[10] = "b5.3"
    arr[11] = arr[12] = "b5.4"
    arr[13] = arr[14] = "B5.1"
    arr[15] = arr[16] = "B5.2"
    return arr[i]
end

function index_to_cell_name_8(i)
    arr = Array{String}(undef,8)
    arr[1] = arr[2] = "B4.1"
    arr[3] = arr[4] = "A4.1"
    arr[5] = arr[6] = "b4.2"
    arr[7] = arr[8] = "a4.2"
    return arr[i]
end

function index_to_cell_name(i,stage)
    if stage == 32
        return index_to_cell_name_32(i)
    elseif stage == 16
        return index_to_cell_name_16(i)
    elseif stage == 8
        return index_to_cell_name_8(i)
    elseif stage == 64
        return index_to_cell_name_64(i)
    else
        error()
    end
end

function index_to_cell_name(i,GE::GeneExpressionStage)
    return index_to_cell_name(i,GE.stage)
end

function index_to_cell_name(i,PC_GE::PCA_GeneExpression)
    return index_to_cell_name(i,PC_GE.GE.stage)
end

function get_contact_graph(stage;all_modes=true)
    if stage ==8
        Ψ1 = [[1 2], [3 4],  [5 6], [7 8]]
        Ψ2 = [[1 3; 2 4], [1 5; 2 6], [3 5; 4 6], [3 7; 4 8], [5 7; 6 8]]
        Ψ_sis = [1 5; 2 6; 1 5; 3 7; 4 8]
    elseif stage == 16
        Ψ1 = [[1 2], [5 6], [7 8],[11 12], [13 14],[15 16]]
        Ψ2 = [[1 3; 2 4], [1 5; 2 6], [1 13; 2 14], [3 5; 4 6], [3 7; 4 8], [3 9; 4 10], [3 13; 4 14], [5 7; 6 8], [7 9; 8 10], [7 11; 8 12], [9 11; 10 12], [9 13; 10 14], [11 13; 12 14], [11 15; 12 16], [13 15; 14 16]]
        Ψ_sis = [1 3; 2 4; 5 7; 6 8; 9 11; 10 12; 13 15; 14 16]
    elseif stage == 32
        Ψ1 = [[1 2], [5 6], [9 10], [11 12], [15 16], [21 22], [23 24], [31 32], [29 30]]
        Ψ2 = [[1 3; 2 4], [1 5; 2 6], [1 29; 2 30], [3 7; 4 8], [3 17; 4 18], [3 15; 4 16], [3 25; 4 26], [3 29; 4 30], [5 7; 6 8], [5 9; 6 10], [5 11; 6 12], [7 9; 8 10], [7 17; 8 18], [7 13; 8 14], [9 11; 10 12], [9 13; 10 14], [13 17; 14 18], [17 19; 18 20], [17 25; 18 26], [11 13; 12 14], [11 15; 12 16], [13 15; 14 16], [13 19; 14 20], [15 19; 16 20], [15 21; 16 22], [19 21; 20 22], [19 25; 20 26], [19 23; 20 24], [21 23; 22 24], [25 27; 26 28], [25 29; 26 30], [23 27; 24 28], [23 31; 24 32], [27 31; 28 32], [27 29; 28 30], [29 31; 30 32]]
        Ψ_sis = [1 5; 2 6; 3 7; 4 8; 9 11; 10 12; 13 15; 14 16; 17 19; 18 20; 21 23; 22 24; 25 29; 26 30; 27 31; 28 32]
    elseif stage == 64
        Ψ1 = [[1 2], [7 8], [9 10], [13 14], [17 18], [35 36], [41 42], [43 44], [49 50], [53 54], [57 58], [59 60], [63 64]]
        Ψ2 = [[1 7; 2 8], [1 5; 2 6], [1 13; 2 14], [5 7; 6 8], [5 25; 6 26], [5 13; 6 14], [5 11; 6 12], [7 9; 8 10], [7 29; 8 30], [7 25; 8 26], [7 13; 8 14], [9 49; 10 50], [9 29; 10 30], [11 25; 12 26], [11 27; 12 28], [11 13; 12 14], [11 15; 12 16], [13 17; 14 18], [13 15; 14 16], [15 27; 16 28], [15 21; 16 22], [15 19; 16 20], [17 35; 18 36], [17 21; 18 22], [19 27; 20 28], [19 31; 20 32], [19 23; 20 24], [19 21; 20 22], [21 23; 22 24], [21 35; 22 36], [21 37; 22 38], [23 33; 24 34], [23 39; 24 40], [23 31; 24 32], [23 37; 24 38], [25 29; 26 30], [25 27; 26 28], [27 29; 28 30], [27 31; 28 32], [29 51; 30 52], [29 47; 30 48], [31 47; 32 48], [31 33; 32 34], [33 3; 34 4], [33 45; 34 46], [33 39; 34 40], [35 41; 36 42], [35 37; 36 38], [37 41; 38 42], [37 39; 38 40], [39 43; 40 44], [39 45; 40 46], [39 41; 40 42], [41 43; 42 44], [43 59; 44 60], [43 61; 44 62], [43 45; 44 46], [45 61; 46 62], [45 3; 46 4], [47 55; 48 56], [47 51; 48 52], [47 3; 48 4], [49 53; 50 54], [49 51; 50 52], [51 53; 52 54], [51 55; 52 56], [53 57; 54 58], [55 57; 56 58], [55 63; 56 64], [55 3; 56 4], [57 63; 58 64], [3 61; 4 62], [59 63; 60 64], [59 61; 60 62], [61 63; 62 64]]
        Ψ_sis = [7 9; 8 10; 27 29; 28 30; 1 5; 2 6; 11 25; 12 26; 31 33; 32 34; 19 23; 20 24; 13 15; 14 16; 17 21; 18 22; 49 53; 50 54; 57 63; 58 64; 47 51; 48 52; 55 3; 56 4; 59 61; 60 62; 41 43; 42 44; 39 45; 40 46; 35 37; 36 38]
    else
        error("Stage is not 8, 16, 32, 64")
    end
    Ψu = [Ψ2...]
    if all_modes
        Ψ_glob = hcat([[i,j] for i = 1:stage for j = (i+1):stage]...) |> permutedims
        Ψ_nbh = vcat(Ψ1...,Ψ2...)
        Ψ_sym = [Ψ_glob,Ψ_sis,Ψ_nbh,Ψ1...]
        
    else
        Ψ_sym = [Ψ1...]
    end
    return Ψ_sym, Ψu
end

function get_contact_graph(GE::GeneExpressionStage;all_modes=true)
    return get_contact_graph(GE.stage,all_modes=all_modes)
end
function get_contact_graph(PC_GE::PCA_GeneExpression;all_modes=true)
    return get_contact_graph(PC_GE.GE.stage,all_modes=all_modes)
end


function print_gene_ranking(r_global,stage,ng;all_modes=false)
    Ψ_sym, ΨU = get_contact_graph(stage,all_modes=all_modes)
    n_el = ((ng+1)*ng)>>1
    nu = ((ng)*ng)
    nΨ = length(Ψ_sym)
    nΨu = length(ΨU)

    id_mat = [(i,j) for i = 1:ng,j=1:ng]
    idx_vec = vcat([vcat([Ψ_sym[i] for k = 1:n_el]) for i = 1:nΨ]...,[vcat([ΨU[i] for k = 1:nu]) for i = 1:nΨu]...)
    gene_vec = vcat([vech(id_mat) for i = 1:nΨ]...,[vcat((id_mat)...) for i = 1:nΨu]...)

    rank_string = String[]
    for i = 1:length(r_global)
        top_coord = findfirst(isequal(i),r_global)
        con = idx_vec[top_coord]
        gg = gene_vec[top_coord]
        if length(con) == stage*(stage-1)
            s = string( i, ". PC ", gg[1]," - PC ", gg[2], " embryo-wide")
        elseif length(con) == stage
            s = string( i, ". PC ", gg[1]," - PC ", gg[2]," sister")
        elseif length(con) > 4
            s = string( i, ". PC ", gg[1]," - PC ", gg[2]," neighbor")
        else
            s = string( i, ". ",index_to_cell_name(con[1,1],stage)," PC ", gg[1]," - ", index_to_cell_name(con[1,2],stage)," PC ",gg[2])
        end
        push!(rank_string,s)
    end
    return rank_string
end
