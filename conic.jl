using JuMP
using Gurobi
using CPLEX
using MathProgBase

require("MIQCPGenerator.jl")



function addTower(threeDimImplementation::Function,epsilon::Float64=0.01)

    function implementTower(model::Model, y, ys, t, threeDimCones::Vector{Vector{JuMP.Variable}})

         dim = length(y)
         K = ceil(log(2,dim))
         skb = ceil(log(4,(16/9)*pi^(-2.0)*log(1+epsilon)))

         currentLevelVars = y
         numCurrentLevelVars = dim
         for k in 1:K-1

             numNextLevelVars = convert(Integer,ceil(numCurrentLevelVars/2.0))
             nextLevelVars = Array(Variable,numNextLevelVars)
             for i in 1:floor(numCurrentLevelVars/2.0)
                 nextLevelVars[i] = Variable(model,0,Inf,:Cont)
             end
             for i in 1:floor(numCurrentLevelVars/2.0)
                 threeDimImplementation(threeDimCones,model,currentLevelVars[2*i-1],currentLevelVars[2*i],nextLevelVars[i],convert(Integer,ceil((k+1)/2)-skb))
             end
             if numCurrentLevelVars%2 == 1
                 nextLevelVars[numNextLevelVars] = currentLevelVars[numCurrentLevelVars]
             end
             numCurrentLevelVars = numNextLevelVars
             currentLevelVars = nextLevelVars
         end
         threeDimImplementation(threeDimCones,model,currentLevelVars[1],currentLevelVars[2],t,convert(Integer,ceil((K+1)/2)-skb))


    end

end

function SeparableThreeDimCone(threeDimCones::Vector{Vector{JuMP.Variable}}, model::Model, x,  y,  t, k )
    @defVar(model, xs >= 0)
    @defVar(model, ys >= 0)
    @addConstraint(model, xs + ys <= t)
    @addConstraint(model, x^2 <= xs*t)
    @addConstraint(model, y^2 <= ys*t)
    push!(threeDimCones,[x,y,t,xs,ys])
end

function StandardThreeDimCone(threeDimCones::Vector{Vector{JuMP.Variable}}, model::Model, x,  y,  t, k )
    push!(threeDimCones,[x,y,t])
    @addConstraint(model, x^2 + y^2 <= t^2)
end

function BNGLIThreeDimCone(threeDimCones::Vector{Vector{JuMP.Variable}}, model::Model, x,  y,  t, k  )
    push!(threeDimCones,[x,y,t])
    alpha = Array(Variable,k+1)
    beta = Array(Variable,k+1)
    alpha[1] = x
    beta[1] = y
    for j in 1:k
        currentCos = cos(pi/(2^(j-1)))
        currentSin = sin(pi/(2^(j-1)))

        alpha[j+1] = Variable(model,-Inf,Inf,:Cont)
        beta[j+1] = Variable(model,-Inf,Inf,:Cont)

        @addConstraint(model, alpha[j+1] == currentCos*alpha[j] + currentSin*beta[j])
        @addConstraint(model, beta[j+1] >= currentCos*beta[j] - currentSin*alpha[j])
        @addConstraint(model, beta[j+1] >= -currentCos*beta[j] + currentSin*alpha[j])
    end
    @addConstraint(model, t == cos(pi/(2^(k)))*alpha[k+1] + sin(pi/(2^(k)))*beta[k+1])
end

function BNGLISepThreeDimCone(threeDimCones::Vector{Vector{JuMP.Variable}}, model::Model, x,  y,  t, k  )
    @defVar(model, xs >= 0)
    @defVar(model, ys >= 0)
    @addConstraint(model, xs + ys <= t)
    push!(threeDimCones,[x,y,t,xs,ys])
    alpha = Array(Variable,k+1)
    beta = Array(Variable,k+1)
    alpha[1] = x
    beta[1] = y
     for j in 1:k
         currentCos = cos(pi/(2^(j-1)))
         currentSin = sin(pi/(2^(j-1)))

         alpha[j+1] = Variable(model,-Inf,Inf,:Cont)
        beta[j+1] = Variable(model,-Inf,Inf,:Cont)

        @addConstraint(model, alpha[j+1] == currentCos*alpha[j] + currentSin*beta[j])
        @addConstraint(model, beta[j+1] >= currentCos*beta[j] - currentSin*alpha[j])
        @addConstraint(model, beta[j+1] >= -currentCos*beta[j] + currentSin*alpha[j])
    end
    @addConstraint(model, t == cos(pi/(2^(k)))*alpha[k+1] + sin(pi/(2^(k)))*beta[k+1])
end

function implementSOCPSeparableBase(model::Model, y, ys, t, threeDimCones::Vector{Vector{JuMP.Variable}})
    dim=length(y)
    @addConstraint(model, sum{ ys[i], i = 1:dim} <= t)
end

function implementSOCPSeparable(model::Model, y, ys, t, threeDimCones::Vector{Vector{JuMP.Variable}})
    dim=length(y)
    @addConstraint(model, sum{ ys[i], i = 1:dim} <= t)
    for i in 1:dim
        @addConstraint(model,y[i]^2<=ys[i]*t)
    end
end

function implementSOCPStandard(model::Model, y, ys, t, threeDimCones::Vector{Vector{JuMP.Variable}} )
    dim=length(y)
    @addConstraint(model, sum{y[i]^2, i=1:dim} <= t^2)
end


function buildModel(prob::MISOCPInput,SOCPImplementations,solver=MathProgBase.defaultQPsolver,equality=true; relaxed = false)
    model = Model(solver=solver)

    nx = size(prob.A,2)
    nz = size(prob.B,2)
    m = size(prob.A,1)

    @defVar(model,prob.lx[i]<=x[i=1:nx]<=prob.ux[i])
    z = []
    if relaxed
        @defVar(model,prob.lz[i]<=z[i=1:nz]<=prob.uz[i])
    else
        @defVar(model,prob.lz[i]<=z[i=1:nz]<=prob.uz[i],Int)
    end

    @setObjective(model, Min, sum{ prob.c[i]*x[i], i = 1:nx} + sum{ prob.f[i]*z[i], i = 1:nz})

    A = prob.A'
    B = prob.B'
    for i=1:m
            if prob.sense[i] == '='
                   @addConstraint(model, sum{ A.nzval[idx]*x[A.rowval[idx]],
                    idx = A.colptr[i]:(A.colptr[i+1]-1)} +sum{ B.nzval[idx]*z[B.rowval[idx]],
                    idx = B.colptr[i]:(B.colptr[i+1]-1)}  == prob.b[i])
               elseif prob.sense[i] == '>'
                    @addConstraint(model, sum{ A.nzval[idx]*x[A.rowval[idx]],
                    idx = A.colptr[i]:(A.colptr[i+1]-1)} +sum{ B.nzval[idx]*z[B.rowval[idx]],
                    idx = B.colptr[i]:(B.colptr[i+1]-1)}  >= prob.b[i])
            else
            @addConstraint(model, sum{ A.nzval[idx]*x[A.rowval[idx]],
                    idx = A.colptr[i]:(A.colptr[i+1]-1)} +sum{ B.nzval[idx]*z[B.rowval[idx]],
                    idx = B.colptr[i]:(B.colptr[i+1]-1)}  <= prob.b[i])
            end
    end

    nsoc = length(prob.D)

    y = Dict()
    ys = Dict()
    @defVar(model, t[1:nsoc] >= 0)
    threeDimCones = Array(Vector{Vector{JuMP.Variable}},nsoc)
    fill!(threeDimCones,[])
    for k in 1:nsoc


        D = prob.D[k]' # transposed
        E = prob.E[k]' # transposed
        d = prob.d[k]
        p = prob.p[k]
        w = prob.w[k]
        q = prob.q[k]

        dim = size(D,2) # dimension of cone


        @defVar(model,ys[k][1:dim]>=0)

        # y = Dx +Ez - d
        @defVar(model, y[k][1:dim])
        # t = p^Tx +w^Tz - q
        
        for i in 1:dim
            if equality
                @addConstraint(model, y[k][i] == sum{ D.nzval[idx]*x[D.rowval[idx]],
                    idx = D.colptr[i]:(D.colptr[i+1]-1)} + sum{ E.nzval[idx]*z[E.rowval[idx]],
                    idx = E.colptr[i]:(E.colptr[i+1]-1)} - d[i] )
            else
                @addConstraint(model, y[k][i] >= sum{ D.nzval[idx]*x[D.rowval[idx]],
                    idx = D.colptr[i]:(D.colptr[i+1]-1)} + sum{ E.nzval[idx]*z[E.rowval[idx]],
                    idx = E.colptr[i]:(E.colptr[i+1]-1)} - d[i] )
                @addConstraint(model, y[k][i] >= sum{ - D.nzval[idx]*x[D.rowval[idx]],
                    idx = D.colptr[i]:(D.colptr[i+1]-1)} + sum{ - E.nzval[idx]*z[E.rowval[idx]],
                    idx = E.colptr[i]:(E.colptr[i+1]-1)} + d[i] )
            end
        end
        @addConstraint(model, t[k] == sum{ p[j]*x[j], j = 1:nx } + sum{ w[j]*z[j], j = 1:nz } - q)
        for implementation in SOCPImplementations
            implementation(model, y[k], ys[k], t[k], threeDimCones[k])
        end
    end
    return model, x, z, y, ys, t, threeDimCones
end



function checkPrimalSolutionQuality(prob::MISOCPInput,x,z)

    nx = size(prob.A,2)
    nz = size(prob.B,2)
    m = size(prob.A,1)
    
     # Bound Infeasibilities
    boundsinf=0.0
    for i in 1:nx
        boundsinf=max(x[i]-prob.ux[i],boundsinf)
        boundsinf=max(prob.lx[i]-x[i],boundsinf)
    end
    for i in 1:nz
        boundsinf=max(z[i]-prob.uz[i],boundsinf)
        boundsinf=max(prob.lz[i]-z[i],boundsinf)
    end
    

    linearinf = 0.0
    A = prob.A'
    B = prob.B'
    for row in 1:m
        lhs = 0.0
        for idx in A.colptr[row]:(A.colptr[row+1]-1)
            lhs+=A.nzval[idx]*x[A.rowval[idx]]
        end
        for idx in B.colptr[row]:(B.colptr[row+1]-1)
            lhs+=B.nzval[idx]*z[B.rowval[idx]]
        end
        if prob.sense[row] == '='
            linearinf = max(abs(lhs-prob.b[row]),linearinf)
        elseif prob.sense[row] == '>'
            linearinf = max(prob.b[row]-lhs,linearinf)
        else
            linearinf = max(lhs-prob.b[row],linearinf)
        end
    end

    quadraticinf = 0.0
    nsoc = length(prob.D)
    for k in 1:nsoc
        D = prob.D[k]' # transposed
        E = prob.E[k]' # transposed
        d = prob.d[k]
        p = prob.p[k]
        w = prob.w[k]
        q = prob.q[k]

        dim = length(d) # dimension of cone
        # y[k,1:dim] = D_k V lambda  - d_k
        # y[k,0] = p_k^T V lambda - q_k
        
        squaresum=0.0

        for i in 1:dim
            temp = -d[i]
            for idx in D.colptr[i]:(D.colptr[i+1]-1)
                temp+=D.nzval[idx]*x[D.rowval[idx]]
            end
            for idx in E.colptr[i]:(E.colptr[i+1]-1)
                temp+=E.nzval[idx]*z[E.rowval[idx]]
            end    
            squaresum += temp^2
        end
        squaresum -= (dot(p,x)+dot(w,z)-q)^2  
        quadraticinf=max(squaresum,quadraticinf)
    end
    return quadraticinf, linearinf, boundsinf


end


function getCorrectedSolution(prob::MISOCPInput, z_lower, z_upper, hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)

    saveIncumbent = true
    for i in 1:length(z_lower)
        setLower(hz[i],z_lower[i])
        setUpper(hz[i],z_upper[i])
        if z_lower[i] < z_upper[i]  - 1e-6
            saveIncumbent = false
        end
    end

    status = solve(hmodel)
    obj = Inf
    if  status == :Optimal && getObjectiveValue(hmodel) < incumbentobj
        obj = getObjectiveValue(hmodel)
        if saveIncumbent
            for i in 1:length(prob.f)
                zstar[i] = getValue(hz[i])
            end
            for i in 1:length(prob.c)
                xstar[i] = getValue(hx[i])
            end
            nsoc = length(prob.D)
            for k in 1:nsoc
                ystar[k] = getValue(hy[k])
                ysstar[k] = getValue(hys[k])
                tstar[k] = getValue(ht[k])
                for i in 1:length(hthreeDimCones[k])
                    threeDimConesstar[k][i]=[getValue(hthreeDimCones[k][i][j]) for j=1:length(hthreeDimCones[k][i])]
                end
            end
        end
    end

    obj
end


function liftedlpsolve(prob::MISOCPInput, masterimplementation, correctionimplementation , mastersolver=MathProgBase.defaultQPsolver, correctionsolver=MathProgBase.defaultQPsolver ,equality=true)
    model, x, z, y, ys, t, threeDimCones = buildModel(prob,masterimplementation,mastersolver,equality)
    hmodel, hx, hz, hy, hys, ht, hthreeDimCones = buildModel(prob,correctionimplementation,correctionsolver,equality; relaxed=true)


    nsoc = length(prob.D)

    newincumbent = false
    incumbentobj = Inf
    xstar = Array(Float64,length(prob.c))
    zstar = Array(Float64,length(prob.f))
    ystar = Array(Any,nsoc)
    ysstar = Array(Any,nsoc)
    tstar = Array(Float64,nsoc)
    threeDimConesstar = []

    function branchcallback(cb)

        feas = cbgetintfeas(cb)
        z_val = getValue(z)

        feasible = true
        for i in 1:length(prob.f)
            if feas[z[i].col] == 1
                feasible = false
                break
            end
        end

        if feasible
            lb = cbgetnodelb(cb)
            ub = cbgetnodelb(cb)
            fixed = true
            firstfractional = 0
            for i in 1:length(prob.f)
                if lb[z[i].col] < ub[z[i].col]  - 1e-6
                    fixed = false
                    firstfractional = i
                    break
                end
            end
            if !fixed
                obj = getCorrectedSolution(prob, [ lb[z[i].col]  for i=1:length(prob.f) ], [ ub[z[i].col]  for i=1:length(prob.f) ], hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)
                if obj < incumbentobj
                    addBranch(cb, z[firstfractional] >= floor(z_val[firstfractional]) + 1, obj)
                    addBranch(cb, z[firstfractional] <= floor(z_val[firstfractional]) , obj)
                end
            end
        end
    end

    function heuristic(cb)

        if newincumbent
            for i in 1:length(prob.c)
                setSolutionValue!(cb, x[i], xstar[i])
            end
            for i in 1:length(prob.f)
                setSolutionValue!(cb, z[i], zstar[i])
            end
            for k in 1:nsoc
                for i in 1:length(ystar[k])
                    setSolutionValue!(cb, y[k][i], ystar[k][i])
                end
                setSolutionValue!(cb, t[k], tstar[k])
            end

            candidateobj = 0
            for i in 1:length(prob.c)
                 candidateobj += prob.c[i]*xstar[i]
            end
            for i in 1:length(prob.f)
                 candidateobj += prob.f[i]*zstar[i]
            end
            addSolution(cb)
            newincumbent = false
        end

    end

    function inccallback(cbdata)
        candidateobj = 0.0
        for i in 1:length(prob.c)
             candidateobj += prob.c[i]*getValue(x[i])
         end
        for i in 1:length(prob.f)
             candidateobj += prob.f[i]*getValue(z[i])
        end

        reject = false
        for k in 1:nsoc
            normy = 0
            for i in 1:length(y[k])
                normy += getValue(y[k][i])^2
            end
            if sqrt(normy) > 1e-6 + getValue(t[k])
                reject = true
                break
            end

        end


        if reject
            obj = getCorrectedSolution(prob, getValue(z), getValue(z), hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)
            if obj < incumbentobj
                incumbentobj = obj
                newincumbent = true
            end
            rejectIncumbent(cbdata)
        else
            acceptIncumbent(cbdata)
        end
    end



    setIncumbentCallback(model, inccallback)
    setBranchCallback(model, branchcallback)
    setHeuristicCallback(model, heuristic)


    return model, x, z, y, ys, t
end

function lazycutsolve(prob::MISOCPInput, masterimplementation, correctionimplementation , mastersolver=MathProgBase.defaultQPsolver, correctionsolver=MathProgBase.defaultQPsolver ,equality=true;  separable = 0)

    model, x, z, y, ys, t, threeDimCones = buildModel(prob,masterimplementation,mastersolver,equality)
    hmodel, hx, hz, hy, hys, ht, hthreeDimCones = buildModel(prob,correctionimplementation,correctionsolver,equality; relaxed=true)


    nsoc = length(prob.D)

    newincumbent = false
    incumbentobj = Inf
    xstar = Array(Float64,length(prob.c))
    zstar = Array(Float64,length(prob.f))
    ystar = Array(Any,nsoc)
    ysstar = Array(Any,nsoc)
    tstar = Array(Float64,nsoc)
    threeDimConesstar = []

    function towerSeparator(cb)
        separated = false
        for k in 1:nsoc
            for i in 1:length(threeDimCones[k])
                y_val = [getValue(threeDimCones[k][i][1]),getValue(threeDimCones[k][i][2])]
                normy = norm(y_val)
                if normy >  1e-6 + getValue(threeDimCones[k][i][3])
                    y_val /= normy
                     @addLazyConstraint(cb, dot(y_val, [threeDimCones[k][i][1],threeDimCones[k][i][2]]) <= threeDimCones[k][i][3] )
                    separated = true
                end
            end
        end

        if separated && dot(prob.c,getValue(x))+dot(prob.f,getValue(z)) < incumbentobj
            z_val = getValue(z)

            obj = getCorrectedSolution(prob, z_val,z_val, hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)

            if obj < incumbentobj
                incumbentobj = obj
                newincumbent = true
            end
        end
    end

    function towerSepSeparator(cb)
        separated = false
        for k in 1:nsoc
            for i in 1:length(threeDimCones[k])
                y_val = [getValue(threeDimCones[k][i][1]),getValue(threeDimCones[k][i][2])]
                ys_val = [getValue(threeDimCones[k][i][4]),getValue(threeDimCones[k][i][5])]
                t_val = getValue(threeDimCones[k][i][3])

                for i in 1:2
                    #a = sqrt((2*y_val[i])^2+(ys_val[i]-t_val)^2)
                    #if a > 1e-6 + t_val +ys_val[i]
                    if y_val[i]^2 > 1e-6 + t_val*ys_val[i]
                        #@addLazyConstraint(cb,  4*(y_val[i]/a) * y[k][i] + ((ys_val[i] - t_val)/a)*(ys[k][i]-t[k])<= t[k] +ys[k][i])
                        a = y_val[i]/t_val
                        @addLazyConstraint(cb,  2*a*y[k][i] - a^2*t[k] <= ys[k][i] )
                        separated = true
                    end
                end
            end
        end

        if separated && dot(prob.c,getValue(x))+dot(prob.f,getValue(z)) < incumbentobj
            z_val = getValue(z)

            obj = getCorrectedSolution(prob, z_val,z_val, hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)

            if obj < incumbentobj
                incumbentobj = obj
                newincumbent = true
            end
        end
    end

    function sepSeparator(cb)
        separated = false
        for k in 1:nsoc
            y_val = getValue(y[k])[:]
            ys_val = getValue(ys[k])[:]
            t_val = getValue(t[k])
            dim = length(prob.d[k])

            for i in 1:dim
                #a = sqrt((2*y_val[i])^2+(ys_val[i]-t_val)^2)
                #if a > 1e-6 + t_val +ys_val[i]
                if y_val[i]^2 > 1e-6 + t_val*ys_val[i]
                    #@addLazyConstraint(cb,  4*(y_val[i]/a) * y[k][i] + ((ys_val[i] - t_val)/a)*(ys[k][i]-t[k])<= t[k] +ys[k][i])
                    a = y_val[i]/t_val
                    @addLazyConstraint(cb,  2*a*y[k][i] - a^2*t[k] <= ys[k][i] )
                    separated = true
                end
            end
        end

        if separated && dot(prob.c,getValue(x))+dot(prob.f,getValue(z)) < incumbentobj
            z_val = getValue(z)

            obj = getCorrectedSolution(prob, z_val,z_val, hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)

            if obj < incumbentobj
                incumbentobj = obj
                newincumbent = true
            end
        end
    end

    function separator(cb)

        separated = false
        for k in 1:nsoc
            y_val = getValue(y[k])[:]
            normy = norm(y_val)

            if normy > 1e-6 + getValue(t[k])
                y_val /= normy
                 @addLazyConstraint(cb, dot(y_val, y[k]) <= t[k] )
                separated = true
            end
        end

        if separated && dot(prob.c,getValue(x))+dot(prob.f,getValue(z)) < incumbentobj
            z_val = getValue(z)

            obj = getCorrectedSolution(prob, z_val,z_val, hz, hx, hy, hys,  ht, hmodel, hthreeDimCones, xstar, zstar, ystar, ysstar, tstar, threeDimConesstar, incumbentobj)

            if obj < incumbentobj
                incumbentobj = obj
                newincumbent = true
            end
        end
    end

    function heuristic(cb)
        if newincumbent
            for i in 1:length(prob.c)
                setSolutionValue!(cb, x[i], xstar[i])
            end
            for i in 1:length(prob.f)
                setSolutionValue!(cb, z[i], zstar[i])
            end

            for k in 1:nsoc
                for i in 1:length(ystar[k])
                    setSolutionValue!(cb, y[k][i], ystar[k][i])
                    if separable == 1
                        setSolutionValue!(cb, ys[k][i], ysstar[k][i])
                    elseif separable == 2
                        for k in 1:nsoc
                            for i in 1:length(hthreeDimCones[k])
                                for j in 1:length(hthreeDimCones[k][i])
                                    setSolutionValue!(cb, threeDimCones[k][i][j], threeDimConesstar[k][i][j])
                                end
                            end
                        end
                    end
                end
                setSolutionValue!(cb, t[k], tstar[k])
            end
            candidateobj = 0
            for i in 1:length(prob.c)
                 candidateobj += prob.c[i]*xstar[i]
             end
            for i in 1:length(prob.f)
                 candidateobj += prob.f[i]*zstar[i]
            end
            addSolution(cb)
            newincumbent = false
        end

    end



    if separable == 0
        setLazyCallback(model, separator)
    elseif separable == 1
        setLazyCallback(model, sepSeparator)
    else
        threeDimConesstar = Array(Vector{Vector{Float64}},nsoc)
        for k in 1:nsoc
            threeDimConesstar[k]=Array(Vector{Float64},length(hthreeDimCones[k]))
        end
        setLazyCallback(model, towerSeparator)
    end

    setHeuristicCallback(model, heuristic)



    return model, x, z, y, ys, t
end

function solveInstance(misocp,results,name,solver;cplexbasesolver=true)

    m, x, z, y, ys, t = solver()
    tic();
    solve(m)
    solvetime = toc();
    quadraticinf, linearinf, boundsinf = checkPrimalSolutionQuality(misocp,getValue(x),getValue(z))
    nodes=0
    bestbound = 0
    if cplexbasesolver
        nodes=CPLEX.getnodecnt(m.internalModel)
        bestbound = CPLEX.getobjbound(m.internalModel)
    else
        nodes=Gurobi.get_node_count(getrawsolver(m.internalModel))
        bestbound = Gurobi.getobjbound(m.internalModel)
    end
    push!(results,[name, getObjectiveValue(m),solvetime,nodes,quadraticinf, linearinf, boundsinf,bestbound,100.0*abs(bestbound-getObjectiveValue(m))/getObjectiveValue(m)])
end


