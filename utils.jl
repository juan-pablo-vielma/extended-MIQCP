using MathProgBase.MathProgSolverInterface
using Clp
using Gurobi

# min  c^Tx
# s.t. Ax sense b
#      l <= x <= u
type LPData
    c::Vector{Float64} # objective vector 
    A::SparseMatrixCSC{Float64,Int64} # constraint matrix
    sense::Vector{Char} # '<', '>', '='
    b::Vector{Float64}
    l::Vector{Float64}
    u::Vector{Float64}
    vartype::Vector{Char} # 'C', 'I'
end

function LPDataFromMPS(mpsfile::String) 

    println(mpsfile)
    m = model(GurobiSolver())
    loadproblem!(m,convert(ASCIIString,mpsfile))

    c = getobj(m)
    A = getconstrmatrix(m)
    nrow,ncol = size(A)
    xlb = getvarLB(m)
    xub = getvarUB(m)
    l = getconstrLB(m)
    u = getconstrUB(m)
    vartype = getvartype(m)

    sense = similar(l,Char)
    b = similar(l)

    for i in 1:nrow
        if l[i] == u[i] && l[i] != Inf
            sense[i] = '='
            b[i] = l[i]
        elseif l[i] > -1e20 && u[i] < 1e20
            error("Range constraints not supported")
        elseif l[i] > -1e20
            sense[i] = '>'
            b[i] = l[i]
        elseif u[i] < 1e20
            sense[i] = '<'
            b[i] = u[i]
        else
            error()
        end
    end
 
    return LPData(c,A,sense,b,xlb,xub,vartype)


end


