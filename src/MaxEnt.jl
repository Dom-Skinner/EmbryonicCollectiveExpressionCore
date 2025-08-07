using LinearAlgebra
using Distributions
using Optim
using Arpack

#@assert length(d.μ) .== nc*ng

function adj_mat(E,nc)
    # Create an adjacency matrix from a list of edges
    M = zeros(nc,nc)
    for i = 1:size(E,1)
        M[E[i,1],E[i,2]] = 1
        M[E[i,2],E[i,1]] = 1
    end
    return M
end

function a_adj_mat(E,nc)
    # Create an adjacency matrix from a list of edges
    M = zeros(nc,nc)
    for i = 1:size(E,1)
        M[E[i,1],E[i,2]] = 1
        M[E[i,2],E[i,1]] = -1
    end
    return M
end

function graph_construct(ed)
    g  = SimpleGraph(maximum(ed))
  for i = 1:size(ed,1)
    add_edge!(g, ed[i,1], ed[i,2])
  end
  return g
end

function rand_sample(M::MaxEntModel,N_sample)
    return [reshape(rand(M.d,1),M.nc,M.ng) |> permutedims for i = 1:N_sample]
end

function fit_independent_max_ent(C,μ,GE::GeneExpressionStage)
    type_ = get_type(GE)
    # Assuming cell types, but no cell-cell interactions
    n_t = [sum(type_ .== i) for i = unique(type_)]
    h = sum([kron(μ[t] , indicator(t,type_) ) for t = 1:length(n_t)])
    Jinv = sum([kron(C[t],diagm(indicator(t,type_))) for t = 1:length(n_t)])
    return MaxEntModel(MvNormal(h, Jinv),GE.stage,length(GE.genename),GE.genename,C,[],[],[],[])
end

function fit_independent_max_ent(C,μ,PC_GE::PCA_GeneExpression)
    type_ = get_type(PC_GE)
    # Assuming cell types, but no cell-cell interactions
    n_t = [sum(type_ .== i) for i = unique(type_)]
    h = sum([kron(μ[t] , indicator(t,type_) ) for t = 1:length(n_t)])
    Jinv = sum([kron(C[t],diagm(indicator(t,type_))) for t = 1:length(n_t)])
    return MaxEntModel(MvNormal(h, Jinv),PC_GE.GE.stage,PC_GE.nPC,"PC_".*string.(1:PC_GE.nPC),C,[],[],[],[])
end


function obj_fun_SU(Mod::MaxEntModel,GE::GeneExpressionStage,C_tot,λ,keep_vec)
    @assert size(Mod.M[1],1) == GE.ng # incase GE and Mod don't match
    return obj_fun_SU(vcat(vech.(Mod.M)...,vech.(Mod.S)...,(Mod.U...)...),GE.ng,GE.stage,Mod.Ψ_sym,Mod.Ψu,C_tot,get_type(GE),λ,keep_vec)
end


function obj_fun_SU(x,ng,nc,Ψv,Ψu,C_tot::CorrelationBounds,cluster,λ,keep_vec)
    
    
    n_el = ((ng+1)*ng)>>1
    n_ael = ((ng-1)*ng)>>1
    n_u = ng*ng
    nt = length(unique(cluster))
    nΨ = length(Ψv)
    nΨu = length(Ψu) 
    
    if length(keep_vec) == 0
        xt = copy(x)
    else
        xt = zeros(length(keep_vec))
        xt[keep_vec] .= x
    end
       
    M = [ivech(xt[(1+n_el*(i-1)):n_el*i]) for i = 1:nt]
    S = [ivech(xt[ (nt*n_el+1 + n_el*(i-1)):(nt+i)*n_el]) for i = 1:nΨ]
    U = [reshape(xt[ (nt+nΨ)*n_el + (i-1)*n_u .+ (1:n_u)],ng,ng) for i = 1:nΨu]


    J =  Symmetric(sum([kron(M[t],diagm(indicator(t,cluster))) for t = 1:nt]) )
    if nΨ > 0 
        J += Symmetric(sum([kron(S[i],adj_mat(Ψv[i],nc)) for i = 1:nΨ]))
    end
    if nΨu > 0
        J += Symmetric(sum([kron(0.5*U[i] + 0.5*U[i]',adj_mat(Ψu[i],nc)) for i = 1:nΨu]) + sum([kron(0.5*U[i] - 0.5*U[i]',a_adj_mat(Ψu[i],nc)) for i = 1:nΨu]))
    end
    
        
    if !isposdef(J)
        return Inf, zeros(size(x))
    end
    Jinv = inv(J)

    @inline idx_set = (i,q) -> nc*(i-1) + q

    
    
    L = λ[1]*sum(abs.(xt[ (nt*n_el+1):((nt+nΨ)*n_el)])) + 
        λ[2]*sum(abs.(xt[ ((nt+nΨ)*n_el + 1):end])) 
    Y = zeros(size(J))

    for t in 1:nt
        cluster_idx = findall(t.==cluster)
        for i = 1:ng,j=1:ng
            
            a = mean([Jinv[idx_set(i,r),idx_set(j,r)] for r in cluster_idx])
            ωt_ij = (a < C_tot.C_low[t][i,j]) ? (C_tot.C_low[t][i,j] - a) : 0
            ωt_ij += (C_tot.C_high[t][i,j] < a) ? (C_tot.C_high[t][i,j] - a) : 0
            for r in cluster_idx
                Y[idx_set(i,r),idx_set(j,r)] -= 2/length(cluster_idx) *ωt_ij
            end
            L = L + ωt_ij^2
        end
    end

    
    for t1 = 1:nt 
        clst_id1 = findall(cluster.==t1)
        nt1 = length(clst_id1)
        for t2 = t1:nt
            clst_id2 = findall(cluster.==t2)
            nt2 = (t1 == t2) ? nt1 - 1 : length(clst_id2) 
            t_id = triangle_inv(t1,t2,nt)
            for i=1:ng,j=1:ng
                if t1 == t2
                    bb = mean([ Jinv[idx_set(i,a),idx_set(j,b)] for  a in clst_id1 for b in setdiff(clst_id1,a)])
                else
                    bb = mean([ Jinv[idx_set(i,a),idx_set(j,b)] for a in clst_id1, b in clst_id2])
                end
                μ_ij = (bb < C_tot.C_cross_low[t_id][i,j]) ? (C_tot.C_cross_low[t_id][i,j] - bb) : 0
                μ_ij += (C_tot.C_cross_high[t_id][i,j] < bb) ? (C_tot.C_cross_high[t_id][i,j] - bb) : 0
                
                for p in clst_id1
                    for q in ((t1 == t2) ? setdiff(clst_id1,p) : clst_id2)
                        Y[idx_set(i,p),idx_set(j,q)] -= 2 / nt1 / nt2 * μ_ij
                    end
                end
                L = L + μ_ij^2
            end
        end
    end
    
    # implementing the following, slightly optimized
    #Y = Jinv*Y*Jinv
    Y_tmp = zeros(size(J))
    mul!(Y_tmp,Y,Jinv)
    mul!(Y,Jinv,Y_tmp)
    
    
    dLdM = [zeros(ng,ng) for i = 1:nt]  
    for t in 1:nt
        for i = 1:ng,j=1:ng
            s = -sum(Y[idx_set(i,p),idx_set(j,p)] for p in findall(t.==cluster))
            dLdM[t][i,j] += s
            if i != j
                dLdM[t][j,i] += s
            end
        end
    end
    dLdM = vcat(vech.(dLdM)...)
     
    #dLdM = zeros(nt*n_el)


    dLdS = [zeros(ng,ng) for i = 1:nΨ]
    for (k,Ψ_) in enumerate(Ψv)
        for i = 1:ng,j=1:ng
            s = -sum(Y[idx_set(i,Ψ_[idx,1]),idx_set(j,Ψ_[idx,2])] + Y[idx_set(i,Ψ_[idx,2]),idx_set(j,Ψ_[idx,1])]  for idx in 1:size(Ψ_,1))
            dLdS[k][i,j] += s
            if i != j
                dLdS[k][j,i] += s
            end
        end
    end
    dLdS = vcat(vech.(dLdS)...)

    dLdU = [zeros(ng,ng) for i = 1:nΨu]
    for (k,Ω) in enumerate(Ψu)
        for i = 1:ng,j=1:ng
            s = -0.5*sum(Y[idx_set(i,Ω[idx,1]),idx_set(j,Ω[idx,2])]  for idx in 1:size(Ω,1)) - 
                0.5*sum(Y[idx_set(j,Ω[idx,2]),idx_set(i,Ω[idx,1])]  for idx in 1:size(Ω,1))
            dLdU[k][i,j] += 2*s
        end
    end
    dLdU = vcat((dLdU...)...)

    ∇L = vcat(dLdM,dLdS,dLdU)

    ∇L[ (nt*n_el+1):((nt+nΨ)*n_el)] .+=  λ[1]*sign.(xt[ (nt*n_el+1):((nt+nΨ)*n_el)])
    ∇L[  ((nt+nΨ)*n_el + 1):end] .+=  λ[2]*sign.(xt[ ((nt+nΨ)*n_el + 1):end])
   

    if length(keep_vec) > 0
        return L, ∇L[keep_vec]
    else
        return L,∇L
    end
end

function find_maxent_params_SU(C_tot::CorrelationBounds,Ψv,Ψu,nc,ng,cluster;λ=[0.,0.],keep_vec=[])
    
    # start with initial guesses for parameters
    nΨ = length(Ψv)
    nΨu = length(Ψu)
    M0 = inv.(C_tot.C)
    S0 = [zeros(ng,ng) for i = 1:nΨ]
    U0 = [zeros(ng,ng) for i = 1:nΨu]

    x0 = vcat(vech.(M0)...,vech.(S0)...,zeros(length(U0)*ng^2)) # vcat U0 manually due to weird StackOverflow error

    if length(keep_vec) > 0
        x0 = x0[keep_vec]
    end

    
    function fg!(F,G,x)
        # do common computations here
        # ...
       
        L,∇L = obj_fun_SU(x,ng,nc,Ψv,Ψu,C_tot,cluster,λ,keep_vec)
       
        if G != nothing
          # code to compute gradient here
          # writing the result to the vector G
          
          G .= ∇L
        end
        if F != nothing
          # value = ... code to compute objective function
          return L
        end
    end
    
    res = optimize(Optim.only_fg!(fg!),
            x0, Optim.LBFGS(),Optim.Options(outer_iterations = 2500, iterations=100_000,show_trace=false))
    println("min is ", Optim.minimum(res))
    if length(keep_vec) > 0
        xs = zeros(length(keep_vec))
        xs[keep_vec] = Optim.minimizer(res)

    else
        xs = Optim.minimizer(res)
    end
    
    n_el = ((ng+1)*ng)>>1
    n_ael = ((ng-1)*ng)>>1
    n_u = ng*ng
    nt = length(C_tot.C)
    Mfinal = [ivech(xs[(1+n_el*(i-1)):n_el*i]) for i = 1:nt]
    Sfinal = [ivech(xs[ (nt*n_el+1 + n_el*(i-1)):(nt+i)*n_el]) for i = 1:nΨ]
    Ufinal = [reshape(xs[ (nt+nΨ)*n_el + (i-1)*n_u .+ (1:n_u)],ng,ng) for i = 1:nΨu]
    return Mfinal,Sfinal,Ufinal
end

function params_to_model(M::Vector{T},S::Vector,U::Vector,Ψv,Ψu,GE::GeneExpressionStage) where T
    cluster = get_type(GE)
    nc = GE.stage


    Jf =  Symmetric(sum([kron(M[t],diagm(indicator(t,cluster))) for t = 1:length(M)]) )
    if length(Ψv) > 0 
        Jf += Symmetric(sum([kron(S[i],adj_mat(Ψv[i],nc)) for i = 1:length(Ψv)]))
    end
    if length(Ψu) > 0
        Jf+= Symmetric(sum([kron(0.5*U[i] + 0.5*U[i]',adj_mat(Ψu[i],nc)) for i = 1:length(Ψu)]) + sum([kron(0.5*U[i] - 0.5*U[i]',a_adj_mat(Ψu[i],nc)) for i = 1:length(Ψu)]))
    end

    if GE.ng == size(M[1],2)
        return MaxEntModel(MvNormal(zeros(size(Jf,1)),inv(Jf)),GE.stage,GE.ng,GE.genename,M,S,U,Ψv,Ψu)
    else # From a PCA, so take ng from m
        return MaxEntModel(MvNormal(zeros(size(Jf,1)),inv(Jf)),GE.stage,size(M[1],2),"PC_".*string.(1:size(M[1],2)),M,S,U,Ψv,Ψu)
    end
end



function exact_mean_corr(M::MaxEntModel,cluster)
    # Given a max ent model, compute the exact cross correlations between cell types.
    ng = M.ng
    nc = M.nc
    nt = length(unique(cluster))

    C_t_cross = [zeros(ng,ng) for t1 = 1:nt for t2 = t1:nt]
    
    Jinv = Matrix(M.d.Σ)
    idx_set = (i,q) -> nc*(i-1) + q

    for t1 = 1:nt 
        for t2 = t1:nt
            t_id = triangle_inv(t1,t2,nt)
            for l=1:ng,m=1:ng
                if t1 == t2
                    clst_id = findall(cluster.==t1)
                    C_t_cross[t_id][l,m] = mean([ Jinv[idx_set(l,a),idx_set(m,b)] for  a in clst_id for b in setdiff(clst_id,a)])
                else
                    C_t_cross[t_id][l,m] = mean( [Jinv[idx_set(l,a), idx_set(m,b)] for a in findall(cluster.==t1), b in findall(cluster.==t2)])
                end
            end
        end
    end
    return C_t_cross
end

function exact_within_corr(M::MaxEntModel,cluster)
    # Given a max ent model, compute the exact correlations within cell types.
    ng = M.ng
    nc = M.nc
    nt = length(unique(cluster))

    C_t = [zeros(ng,ng) for t1 = 1:nt]
    
    Jinv = Matrix(M.d.Σ)
    idx_set = (i,q) -> nc*(i-1) + q

    for t = 1:nt 
        for l=1:ng,m=1:ng
            C_t[t][l,m] = mean( [Jinv[idx_set(l,a), idx_set(m,a)] for a in findall(cluster.==t)])
        end
    end
    return C_t
end
