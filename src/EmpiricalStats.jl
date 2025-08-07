

nanmean(x) = mean(filter(!isnan,x))


function empirical_mean_corr(GE::GeneExpressionStage;cell_types=true,N_trial=1_000)
    type_ = get_type(GE)
    if !cell_types
        type_ = ones(Int64,size(type_))
    end

    Z = vcat([arrange_GE_data(GE) for idx = 1:N_trial]...)
    C_t,C_t_cross,μ_t = empirical_mean_corr(Z,type_)

    return C_t, C_t_cross, μ_t
end

function empirical_mean_corr(PC_GE::PCA_GeneExpression,N_trial=1_000)
    type_ = get_type(PC_GE.GE)

    Z = vcat([arrange_GE_data(PC_GE) for idx = 1:N_trial]...)
    C_t,C_t_cross,μ_t = empirical_mean_corr(Z,type_)

    return C_t, C_t_cross, μ_t
end


function bootstrap_correlation(GE;N_trial=1_000,return_raw=false,Δ=0.1)
    type_ = get_type(GE)


    C_tot = []
    C_cross_tot = []
    M = pmap(x->empirical_mean_corr(arrange_GE_data(GE,bootstrap=true),type_),1:N_trial)
    for idx = 1:N_trial
        C_t,C_t_cross,_ = M[idx]
        push!(C_tot,vcat(C_t...))
        push!(C_cross_tot,vcat(C_t_cross...))
    end

#=
    C_tot = []
    C_cross_tot = []
    for idx = 1:N_trial
        C_t,C_t_cross,_ = empirical_mean_corr(arrange_GE_data(GE,bootstrap=true),type_)
        push!(C_tot,vcat(C_t...))
        push!(C_cross_tot,vcat(C_t_cross...))
    end
=#
    nt = length(unique(type_))
    ng = find_ng(GE)
    nest_arr = x -> [x[ (1 + (i-1)*ng): i*ng,:] for i = 1:nt]
    nest_arr_cross = x -> [x[ (1 + (i-1)*ng): i*ng,:] for i = 1:((nt*(nt+1))>>1)]
    
    if return_raw
        return [nest_arr(c) for c in C_tot], [nest_arr_cross(c) for c in C_cross_tot]
    end
    
    C_low = copy(C_tot[1])
    C_med = copy(C_tot[1])
    C_high = copy(C_tot[1])
    for i = 1:size(C_low,1),j=1:size(C_low,2)
        C_low[i,j] = quantile(filter(!isnan,[C_tot[k][i,j] for k = 1:length(C_tot)]),Δ/2)
        C_med[i,j] = quantile(filter(!isnan,[C_tot[k][i,j] for k = 1:length(C_tot)]),1/2)
        C_high[i,j] = quantile(filter(!isnan,[C_tot[k][i,j] for k = 1:length(C_tot)]),1-Δ/2)
    end
    
    C_cross_low = copy(C_cross_tot[1])
    C_cross_med = copy(C_cross_tot[1])
    C_cross_high = copy(C_cross_tot[1])
    for i = 1:size(C_cross_low,1),j=1:size(C_cross_low,2)
        C_cross_low[i,j] = quantile(filter(!isnan,[C_cross_tot[k][i,j] for k = 1:length(C_cross_tot)]),Δ/2)
        C_cross_med[i,j] = quantile(filter(!isnan,[C_cross_tot[k][i,j] for k = 1:length(C_cross_tot)]),1/2)
        C_cross_high[i,j] = quantile(filter(!isnan,[C_cross_tot[k][i,j] for k = 1:length(C_cross_tot)]),1-Δ/2)
    end
    C_tot = CorrelationBounds(nest_arr(C_med),nest_arr_cross(C_cross_med),nest_arr(C_low),nest_arr(C_high),nest_arr_cross(C_cross_low),nest_arr_cross(C_cross_high))
    return C_tot

end


function empirical_mean_corr(Z::AbstractVector{T},cluster) where T
      
    ng = size(Z[1],1)
    nt = length(unique(cluster))

    C_t = [zeros(ng,ng) for idx = 1:nt]
    C_t_cross = [zeros(ng,ng) for t1 = 1:nt for t2 = t1:nt]
    
    μ_t = [zeros(ng) for idx = 1:nt]
    for t in 1:nt
        for j=1:ng
            μ_t[t][j] = 0*nanmean([ Z[i][j,a] for i = 1:length(Z), a in findall(cluster.==t)])
        end
    end


    for t = 1:nt
        for l=1:ng,m=1:ng
            C_t[t][l,m] = nanmean([ (Z[i][l,a] - μ_t[t][l]) *(Z[i][m,a] - μ_t[t][m])  for i = 1:length(Z), a in findall(cluster.==t)])
        end

        ##### check for bad stats
        not_nan_cnt = sum(.!isnan.([Z[i][1,a] *Z[i][1,a]  for i = 1:length(Z), a in findall(cluster.==t)]))
        if not_nan_cnt <= ng
                #println(length(Z))
                #println(not_nan_cnt)
		println("Warning: Bad statistics in C_t computation")
		if !any(isnan.(C_t[t]))
		  # if C_t[t] contains NaN's that means there were no cells so just pass a NaN matrix
                  # if C_t[t] exists but has poor stats, make it pos definite
		  C_t[t] .= pos_def_cone(C_t[t]) .+ (1e-8)*I(ng)
        	end
        end
        #@assert isposdef(C_t[t])

    end

    
    for t1 = 1:nt 
        for t2 = t1:nt
            t_id = triangle_inv(t1,t2,nt)
            for l=1:ng,m=1:ng
                if t1 == t2
                    clst_id = findall(cluster.==t1)
                    C_t_cross[t_id][l,m] = nanmean([ (Z[i][l,a] -μ_t[t1][l] ) *(Z[i][m,b] -μ_t[t2][m] ) for i = 1:length(Z), a in clst_id for b in setdiff(clst_id,a)])
                else
                    C_t_cross[t_id][l,m] = nanmean([ (Z[i][l,a] -μ_t[t1][l] ) *(Z[i][m,b] -μ_t[t2][m] ) for i = 1:length(Z), a in findall(cluster.==t1), b in findall(cluster.==t2)])
                end
            end
 
        end
    end

    C_t_cross  = [C_t_cross[triangle_inv(t1,t2,nt)]  for t1 = 1:nt for t2 = t1:nt]
    return C_t,C_t_cross,μ_t
end




function bootstrap_4_correlation(GE;N_trial=1_000,Δ=0.1)
    type_ = get_type(GE)


    M = map(x->empirical_4_point(arrange_GE_data(GE,bootstrap=true),type_),1:N_trial)
    M = hcat(M...)
    
    
    C_low = zeros(size(M,1))
    C_med = zeros(size(M,1))
    C_high = zeros(size(M,1))
    for i = 1:size(C_low,1)
        C_low[i]  = quantile(filter(!isnan,M[i,:]),Δ/2)
        C_med[i]  = quantile(filter(!isnan,M[i,:]),1/2)
        C_high[i] = quantile(filter(!isnan,M[i,:]),1-Δ/2)
    end
    
    
    return C_low, C_med, C_high

end


function bootstrap_3_correlation(GE;N_trial=1_000,Δ=0.1)
    type_ = get_type(GE)


    M = map(x->empirical_3_point(arrange_GE_data(GE,bootstrap=true),type_),1:N_trial)
    M = hcat(M...)
    
    
    C_low = zeros(size(M,1))
    C_med = zeros(size(M,1))
    C_high = zeros(size(M,1))
    for i = 1:size(C_low,1)
        C_low[i]  = quantile(filter(!isnan,M[i,:]),Δ/2)
        C_med[i]  = quantile(filter(!isnan,M[i,:]),1/2)
        C_high[i] = quantile(filter(!isnan,M[i,:]),1-Δ/2)
    end
    
    
    return C_low, C_med, C_high

end
function number_type_3(nt)
    # is there a formula for this number - yes
    # am I going to work it out by hand - no, this will do
    idx = 0
    for t1 = 1:nt 
        for t2 = (t1+1):nt
            for t3 = (t2+1):nt
                idx += 1
            end
        end
    end
    return idx
end

function empirical_3_point(Z::AbstractVector{T},cluster) where T
    nanmean(x) = mean(filter(!isnan,x))
    ng = size(Z[1],1)
    nt = length(unique(cluster))
    n_nt = number_type_3(nt)
    C_t_cross = zeros(n_nt,ng,ng,ng)

    
    idx=0
    for t1 = 1:nt 
        idx1 = findall(cluster.==t1)
        for t2 = (t1+1):nt
            idx2 = findall(cluster.==t2)
            for t3 = (t2+1):nt
                idx3 = findall(cluster.==t3)
                idx += 1
                for l=1:ng,m=1:ng,n=1:ng
                    C_t_cross[idx,l,m,n] = nanmean([ Z[i][l,a] *Z[i][m,b]*Z[i][n,c] for i = 1:length(Z),a in idx1, b in idx2, c in idx3]) 
                end
            end
        end
    end
    return vcat(C_t_cross...)
end


function number_type_4(nt)
    # is there a formula for this number - yes
    # am I going to work it out by hand - no, this will do
    idx = 0
    for t1 = 1:nt 
        for t2 = (t1+1):nt
            for t3 = (t2+1):nt
                for t4 = (t3+1):nt
                    idx += 1
                end
            end
        end
    end
    return idx
end

function empirical_4_point(Z::AbstractVector{T},cluster) where T
    nanmean(x) = mean(filter(!isnan,x))
    ng = size(Z[1],1)
    nt = length(unique(cluster))
    n_nt = number_type_4(nt)
    C_t_cross = zeros(n_nt,ng,ng,ng,ng)

    
    idx=0
    for t1 = 1:nt 
        idx1 = findall(cluster.==t1)
        for t2 = (t1+1):nt
            idx2 = findall(cluster.==t2)
            for t3 = (t2+1):nt
                idx3 = findall(cluster.==t3)
                for t4 = (t3+1):nt
                    idx4 = findall(cluster.==t4)
                    idx += 1
                    for l=1:ng,m=1:ng,n=1:ng,o=1:ng 
                        C_t_cross[idx,l,m,n,o] = nanmean([ Z[i][l,a] *Z[i][m,b]*Z[i][n,c] *Z[i][o,d] for i = 1:length(Z),a in idx1, b in idx2, c in idx3, d in idx4])
                       
                    end
                end
            end
        end
    end
    return vcat(C_t_cross...)
end



function exact_4_point(Mod,cluster)
    
    ng = Mod.ng
    nc = Mod.nc
    nt = length(unique(cluster))
    idx_set = (i,q) -> nc*(i-1) + q
    n_nt = number_type_4(nt)
    C_t_cross = zeros(n_nt,ng,ng,ng,ng)

    
    idx=0
    for t1 = 1:nt 
        idx1 = findall(cluster.==t1)
        for t2 = (t1+1):nt
            idx2 = findall(cluster.==t2)
            for t3 = (t2+1):nt
                idx3 = findall(cluster.==t3)
                for t4 = (t3+1):nt
                    idx4 = findall(cluster.==t4)
                    idx += 1
                    for l=1:ng,m=1:ng,n=1:ng,o=1:ng 
                        C_t_cross[idx,l,m,n,o] = mean([ Mod.d.Σ[idx_set(l,a),idx_set(m,b)]*Mod.d.Σ[idx_set(n,c),idx_set(o,d)] +
                                                        Mod.d.Σ[idx_set(l,a),idx_set(n,c)]*Mod.d.Σ[idx_set(m,b),idx_set(o,d)] +
                                                        Mod.d.Σ[idx_set(l,a),idx_set(o,d)]*Mod.d.Σ[idx_set(m,b),idx_set(n,c)]  for a in idx1, b in idx2, c in idx3, d in idx4])
                    end
                end
            end
        end
    end
    return vcat(C_t_cross...)
end


function Gaussian_4_point(C,cluster)
    ng = size(C.C[1],1)
    nt = length(unique(cluster))
    n_nt = Ascidian.number_type_4(nt)
    C_t_cross = zeros(n_nt,ng,ng,ng,ng)

    
    idx=0
    for t1 = 1:nt 
        for t2 = (t1+1):nt
            for t3 = (t2+1):nt
                for t4 = (t3+1):nt
                    idx += 1
                    for l=1:ng,m=1:ng,n=1:ng,o=1:ng 
                        C_t_cross[idx,l,m,n,o] = C.C_cross[triangle_inv(t1,t2,nt)][l,m]* C.C_cross[triangle_inv(t3,t4,nt)][n,o] + 
                                                    C.C_cross[triangle_inv(t1,t3,nt)][l,n]* C.C_cross[triangle_inv(t2,t4,nt)][m,o] + 
                                                    C.C_cross[triangle_inv(t1,t4,nt)][l,o]* C.C_cross[triangle_inv(t2,t3,nt)][m,n] 
                    end
                end
            end
        end
    end
    return vcat(C_t_cross...)
end