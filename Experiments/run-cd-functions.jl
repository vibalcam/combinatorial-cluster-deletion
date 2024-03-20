"""
Code to run each of the methods independently
"""

using Random

seed = 123456

"""
Run MFP (cluster deletion) for one graph A.
"""
function run_matchflippivot_cd(A,Elist,pivtimes)
    Random.seed!(seed)
    # Get lower bounds
    # tic = time()
    timer = @elapsed bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)
    # timer = time()-tic

    # Round the matching
    Random.seed!(seed)
    # tic = time()
    timer_r = @elapsed r_match, clus, comp, timings, objectives = CD_LP_round_timed(A,Elist,soln_match,pivtimes, true)
    # timer_r = time()-tic
    # println("Upper bound done")
    # @assert(isCDfeasible(A,clus))

    # println("Done checking feasibility for CD ")
    runtimes = [timer, timer_r, timings]
    # timer: time to run MFP
    # timer_r: time to round
    results = [bd, r_match, comp, objectives]
    # bd: lower bound
    # r_match: upper bound

    return clus, runtimes, results
end

"""
Run deterministic MFP (cluster deletion) for one graph A.
"""
function run_det_matchflippivot_cd(A,Elist,pivtimes, degree_based = false)
    # if pivtimes is not one, raise an error
    if pivtimes != 1
        throw(Exception("pivtimes must be 1 for deterministic MFP"))
    end

    # Get lower bounds
    Random.seed!(seed)
    # tic = time()
    timer = @elapsed bd, Elist_match, soln_match = maximal_disjoint_openwedge_rounding(A)
    # timer = time()-tic

    # number of deleted edges (E_W)
    E_1 = sum(soln_match)
    E_0 = length(soln_match) - E_1

    # Round the matching
    Random.seed!(seed)
    # tic = time()
    if degree_based
        timer_r = @elapsed clus = det_pivot_degree(A,Elist,soln_match)
    else
        # clus = cd_deterministic_rounding(A,Elist,soln_match)
        timer_r = @elapsed clus = det_pivot_ratio(A,Elist,soln_match)
    end
    # timer_r = time()-tic
    r_match, M_0, M_h, M_1 = check_cd_obj_errors(Elist,clus, soln_match)
    # println("Upper bound done")

    # checks
    # @assert(isCDfeasible(A,clus))
    @assert(M_h == 0)
    @assert(M_1 + M_0 + M_h == r_match)
    # also check bound on number of mistakes
    println("M_0/E_0: $(M_0/E_0)")
    println("M_1/E_1: $(M_1/E_1)")
    @assert(M_1 >= E_1/2)

    # println("Done checking feasibility for CD")
    runtimes = [timer, timer_r]
    # timer: time to run MFP
    # timer_r: time to round
    results = [bd, r_match, -1]
    # bd: lower bound
    # r_match: upper bound

    return clus, runtimes, results
end

"""
Run min st cut LP (cluster deletion) for one graph A.
"""
function run_stcminst(A, Elist, pivtimes, tl, det_rounding=false, use_min_source = false)
    Random.seed!(seed)
    # tic = time()
    timer = @elapsed bd, soln_stclp, time_st = CDSTC_MINST(A, Elist, use_min_source)
    # time_st: runtime min s-t solver (once defined graph for min s-t)
    # timer = time()-tic

    # number of deleted edges (E_W)
    E_1 = sum(soln_stclp .== 1)
    E_0 = sum(soln_stclp .== 0)
    # to compute number of 1/2 substract number of 1 and 0 from total number of edges to no have to deal with tolerances
    E_h = length(soln_stclp) - E_0 - E_1

    # # check objective
    # stat_stclp, bd_lp, _, _ = CDSTC_Gurobi(A; clusterdeletion = false, LP = true,timelimit = tl,outputflag = false)
    # # Same as CD-LP?
    # print("bd: $bd, bd_lp: $bd_lp, stat_stclp: $stat_stclp")
    # if stat_stclp && bd != bd_lp 
    #     throw(Exception("Unexpected optimal value"))
    # end
    iscdlp = check_cd_constraints(A,Elist,soln_stclp)

    # Round STCminst
    Random.seed!(seed)
    # tic = time()
    if det_rounding
        timer_r = @elapsed clus = det_pivot_degree(A,Elist,soln_stclp)
        r_stclp, M_0, M_h, M_1 = check_cd_obj_errors(Elist,clus, soln_stclp)
        comp = 0
    else
        timer_r = @elapsed r_stclp, clus, comp = CD_LP_round(A,Elist,soln_stclp,pivtimes)

        r_check, M_0, M_h, M_1 = check_cd_obj_errors(Elist,clus, soln_stclp)
        @assert(r_check == r_stclp)
    end
    # timer_r = time()-tic

    # checks
    # @assert(isCDfeasible(A,clus))
    @assert(M_1 + M_0 + M_h == r_stclp)
    # also check bound on number of mistakes
    println("M_0/E_0: $(M_0/E_0)")
    println("M_h/E_h: $(M_h/E_h)")
    println("M_1/E_1: $(M_1/E_1)")
    @assert(M_h >= E_h/2)

    runtimes = [timer, timer_r]
    # timer: time to run min st method
    # timer_r: time to round
    results = [bd, r_stclp, comp, true, iscdlp, time_st]
    # bd: lower bound
    # r_stclp: upper bound

    return clus, runtimes, results
end

"""
Run STC LP (cluster deletion) for one graph A.
"""
function run_stclp(A,Elist,pivtimes,tl,det_rounding=false)
    Random.seed!(seed)
    # tic = time()
    timer = @elapsed stat_stclp, bd, Elist_stc, soln_stclp = CDSTC_Gurobi(A; clusterdeletion = false, LP = true, timelimit = tl,outputflag = false)
    # timer = time()-tic

    # Same as CD-LP?
    iscdlp = check_cd_constraints(A,Elist,soln_stclp)

    # Round STCLP
    Random.seed!(seed)
    # tic = time()
    if det_rounding
        timer_r = @elapsed clus = det_pivot_degree(A,Elist,soln_stclp)
        r_stclp, M_0, M_h, M_1 = check_cd_obj_errors(Elist,clus, soln_stclp)
        comp = 0
    else
        timer_r = @elapsed r_stclp, clus, comp = CD_LP_round(A,Elist,soln_stclp,pivtimes)
    end
    # timer_r = time()-tic

    # if stat_stclp
    #     # If here, then the LP solver converged
    #     # @assert(isCDfeasible(A,clus))
    # end

    runtimes = [timer, timer_r]
    # timer: time to run stc lp
    # timer_r: time to round
    results = [bd, r_stclp, comp, stat_stclp, iscdlp]

    return clus, runtimes, results
end

# # 2-approximation, full cluster deletion LP relaxation
# function run_cdlp(A,Elist,pivtimes,tl)
#     tic = time()
#     stat_cdlp, bd, Elist_cd, soln_cdlp = CDSTC_Gurobi(A; clusterdeletion = true, LP = true,timelimit = tl,outputflag = false)
#     timer = time()-tic
#     @assert(Elist_cd == Elist)

#     tic = time()
#     r_cdlp, clus, comp = CD_LP_round(A,Elist,soln_cdlp,pivtimes)
#     timer_r = time()-tic
#     if stat_cdlp
        # @assert(isCDfeasible(A,clus))
#     end

#     runtimes = [timer, timer_r]
#     results = [bd, r_cdlp, comp, stat_cdlp]

#     return clus, runtimes, results
# end

# # Find optimal solution to the minSTC problem
# function run_stc_opt(A,Elist,tl)
#     tic = time()
#     stat_stc, stc_opt, Elist_stc, soln_stc_opt = CDSTC_Gurobi(A; clusterdeletion = false, LP = false,timelimit = tl)
#     runtime = time()-tic
#     @assert(Elist_stc == Elist)

#     # Check if answer also is same as CD-ILP
#     iscd = check_cd_constraints(A,Elist,soln_stc_opt)
    
#     # Can be used as a 2-approx for CD.
#     r_stc_opt, clus, comp = CD_LP_round(A,Elist,soln_stc_opt,1)

#     results = [stc_opt, r_stc_opt, comp, stat_stc, iscd]

#     return clus, runtime, results
# end

# # Find optimal solution to cluster deletion 
# function run_cd_opt(A,Elist,tl)
#     tic = time()
#     stat_cd, cd_opt, Elist_cd, soln_cd_opt = CDSTC_Gurobi(A; clusterdeletion = true, LP = false,timelimit = tl)
#     runtime = time()-tic
#     @assert(Elist_cd == Elist)

#     r_cd_opt, clus, comp = CD_LP_round(A,Elist,soln_cd_opt,1)
    
#     if stat_cd
#         @assert(round(Int64,r_cd_opt) == round(Int64,cd_opt))
#     end

#     results = [cd_opt, comp, stat_cd]

#     return clus, runtime, results
# end
