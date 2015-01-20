using JuMP


# min  c^Tx +f^Tz 
# s.t. Ax +B z sense b
#      ||D_ix + E_i z- d_i|| <= p_i^Tx + w_i^T z- q_i, i = 1, ..., k
#      lx <= x <= ux
#      lz <= z <= uz
#            z_i integer for all i
# sense is vector of '<', '=', '>'

type MISOCPInput
    c::Vector{Float64}
    f::Vector{Float64}	
    A::SparseMatrixCSC{Float64,Int}
    B::SparseMatrixCSC{Float64,Int}
    b::Vector{Float64}
    sense::Vector{Char}
    D::Vector{SparseMatrixCSC{Float64,Int}}
    E::Vector{SparseMatrixCSC{Float64,Int}}
    d::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
    w::Vector{Vector{Float64}}	
    q::Vector{Float64}
    lx::Vector{Float64} # -Inf means no bound
    ux::Vector{Float64} # Inf means no bound
    lz::Vector{Float64} # -Inf means no bound
    uz::Vector{Float64} # Inf means no bound
end



#################################################
# Shortfall Portfolio
#
# * Indicate Constraints dropped if MIP=false 
#
# min_{x}   - r x
# 
# 1x = 1 				
# ||p1 covHalf x||_2 <= r x - l1 
# ||p2 covHalf x||_2 <= r x - l2	
# x_i <= z_i 				*
# 1z <= k 					*
# x >= 0
# z_i in {0,1} 				*
#
# The investment options include a riskless 
# asset
#
# Note that because of numerical issues with 
# p1*covHalf and p2*covHalf and the precision 
# of .mps file output, the generated problems
# are slightly different to the original .mps
# files. 
#
################################################



function buildShortfallMISOCP(r,covHalf,p1,l1,p2,l2,MIP=true,k=10)

	
	n = length(r)+1
	c = -[r[:]+1;1]

	d = Array(Vector{Float64}, 2)
	d[1] = zeros(n-1)
	d[2] = zeros(n-1)
	p = Array(Vector{Float64}, 2)
	p[1] = [r[:]+1;1]
	p[2] = [r[:]+1;1]
	D = Array(SparseMatrixCSC{Float64,Int}, 2)
	D[1] = sparse([p1*covHalf zeros(n-1)])
	D[2] = sparse([p2*covHalf zeros(n-1)])
	E = Array(SparseMatrixCSC{Float64,Int}, 2)
	w = Array(Vector{Float64}, 2)
	q = Array(Float64, 2)
	q[1] = l1
	q[2] = l2

	lx = zeros(n)
	ux = ones(n)
	if MIP
		f = zeros(n)
		A = [speye(n); sparse(ones(1,n)); spzeros(1,n)]
	 	B = [-speye(n); spzeros(1,n); sparse(ones(1,n))]		       
		b = [zeros(n);1.0; k]
		sense = [['<' for i in [1:n]],'=','<']

		E[1] = zeros(n-1,n)
		w[1] = zeros(n)
		E[2] = zeros(n-1,n)
		w[2] = zeros(n)

		lz = zeros(n)
		uz = ones(n)
	else
		f = Float64[]
		A = sparse([ones(1,n)])	
		B = spzeros(1,0)
		b = [1.0]
		sense = ['=']

		E[1] = zeros(n-1,0)
		w[1] = zeros(0)
		E[2] = zeros(n-1,0)
		w[2] = zeros(0)

		lz = Float64[]
		uz = Float64[]
	end
	misocp = MISOCPInput(c, f, A, B, b, sense, D, E, d, p, w, q, lx, ux, lz, uz)
	return misocp
end

function buildShortfallPor(porfile::String,MIP=true)
	r, covHalf=loadPorFile(porfile)
	return buildShortfallMISOCP(r,covHalf,0.841621,0.9,1.88079,0.7,MIP)
end


#################################################
# Standard Markowitz Portfolio
#
# * Indicate Constraints dropped if MIP=false 
#
# min_{x}   - r x
# 
# 1x = 1 					
# ||covHalf x||_2 <= s 	
# x_i <= z_i 				*
# 1z <= k 					*
# x >= 0
# z_i in {0,1} 				*
#
################################################


function buildMarkowitzMISOCP(r,covHalf,s=0.2,MIP=true,k=10)

	n = length(r)
	c = -r[:]

	d = Array(Vector{Float64}, 1)
	d[1] = zeros(n)
	p = Array(Vector{Float64}, 1)
	p[1] = zeros(n)
	D = Array(SparseMatrixCSC{Float64,Int}, 1)
	D[1] = sparse(covHalf)
	E = Array(SparseMatrixCSC{Float64,Int}, 1)
	w = Array(Vector{Float64}, 1)
	q = Array(Float64, 1)
	q[1] = -s


	lx = zeros(n)
	ux = ones(n)
	if MIP
		f = zeros(n)
		A = [speye(n); sparse(ones(1,n)); spzeros(1,n)]
	 	B = [-speye(n); spzeros(1,n); sparse(ones(1,n))]		       
		b = [zeros(n);1.0; k]
		sense = [['<' for i in [1:n]],'=','<']

		E[1] = zeros(n,n)
		w[1] = zeros(n)

		lz = zeros(n)
		uz = ones(n)
	else
		f = Float64[]
		A = sparse([ones(1,n)])	
		B = spzeros(1,0)
		b = [1.0]
		sense = ['=']

		E[1] = zeros(n,0)
		w[1] = zeros(0)

		lz = Float64[]
		uz = Float64[]
	end
	misocp = MISOCPInput(c, f, A, B, b, sense, D, E, d, p, w, q, lx, ux, lz, uz)
	return misocp
end

function buildMarkowitzPor(porfile::String,MIP=true)
	r, covHalf=loadPorFile(porfile)
	return buildMarkowitzMISOCP(r,covHalf,0.2,MIP)
end

function loadPorFile(porfile::String)
	# Data = readdlm(porfile,' ')
	# n = int(Data[1,1])
	# covHalf= float(Data[3:(2+n),1:n])	
	# r = Data[2,1:n]

	# temp=diag(covHalf*covHalf')
	# @printf("Fraction of assets bellow risk threshold = %f\n", length(find(x->x<=0.2^2,temp))/n)
	# return r[:], covHalf
	file = open(porfile,"r")
	n = int(readline(file))
	r = float(split(readline(file))[1:n])
	covHalf = zeros(n,n)
	for i in 1:n
		covHalf[i,:]=float(split(readline(file))[1:n])
	end
	temp=diag(covHalf*covHalf')
	@printf("Fraction of assets bellow risk threshold = %f\n", length(find(x->x<=0.2^2,temp))/n)
	return r[:], covHalf
end	



#####################################################
# Standard Markowitz Portfolio with Robust Objective
#
# * Indicate Constraints dropped if MIP=false 
#
# min_{x,t}   - t
# 
# 1x = 1 	
# ||covHalf x||_2 <= s 	
# ||3 RHalf x||_2 <= a x - t   				
# x_i <= z_i 				*
# 1z <= k 					*
# x >= 0
# z_i in {0,1} 				*
###################################################

function buildRobustMarkowitzMISOCP(a, RHalf, covHalf,s=0.2,MIP=true,k=10)

	n = length(a)
	c = [zeros(n); -1]

	d = Array(Vector{Float64}, 2)
	d[1] = zeros(n) 
	d[2] = zeros(n) 
	p = Array(Vector{Float64}, 2)
	p[1] = zeros(n+1)
	p[2] = [a[:]; -1]
	D = Array(SparseMatrixCSC{Float64,Int}, 2)
	D[1] = sparse([covHalf zeros(n)])
	D[2] = sparse([3*RHalf zeros(n)]) 
	E = Array(SparseMatrixCSC{Float64,Int}, 2)
	w = Array(Vector{Float64}, 2)
	q = Array(Float64, 2)
	q[1] = -s
	q[2] = 0


	lx = [zeros(n); -Inf]  
	ux = [ones(n); Inf] 
	if MIP
		f = zeros(n) 
		A = [speye(n) spzeros(n,1); sparse([ones(1,n) 0]); spzeros(1,n+1)]
	 	B = [-speye(n); spzeros(1,n); sparse(ones(1,n))]		       
		b = [zeros(n); 1.0; k]
		sense = [['<' for i in [1:n]],'=','<']

		E[1] = zeros(n,n) 
		E[2] = zeros(n,n)
		w[1] = zeros(n)
		w[2] = zeros(n)

		lz = zeros(n)
		uz = ones(n)
	else
		f = Float64[]
		A = sparse([ones(1,n) 0])	
		B = spzeros(1,0)
		b = [1.0]
		sense = ['=']

		E[1] = zeros(n,0)
		E[2] = zeros(n,0)
		w[1] = zeros(0)
		w[2] = zeros(0)

		lz = Float64[]
		uz = Float64[]
	end
	misocp = MISOCPInput(c, f, A, B, b, sense, D, E, d, p, w, q, lx, ux, lz, uz)
	return misocp
end

#####################################################
# Standard Markowitz Portfolio with Robust Objective
# Version for original .por files
#  
# The code used to generate the original .mps files
# from the .por files had some formulation artifacts
# that results in a slightly different problem. 
# This problem is equivalent to the one generated
# by buildRobustMarkowitzMISOCP for the data in the
# .por files, but may yield different alternative
# optimal solutions and/or solution times. 
#
# buildRobustMarkowitzPorMISOCP generates a version
# of the problem that more closely matches the .mps 
# files. 
#
# * Indicate Constraints dropped if MIP=false 
#
# min_{x,t}   - t
# 
# 1x = 1 	
# ||covHalf x||_2 <= s 	
# ||3 covHalf x||_2 <= a x - t   				
# x_i <= z_i 				*
# t <= z_(n+1)
# 1z <= k +1				*
# x >= 0
# t >= 0
# z_i in {0,1} 				*
###################################################


function buildRobustMarkowitzOriginalPorMISOCP(a, RHalf, covHalf,s=0.2,MIP=true,k=10)

	n = length(a)
	c = [zeros(n); -1]

	d = Array(Vector{Float64}, 2)
	d[1] = zeros(n) 
	d[2] = zeros(n) 
	p = Array(Vector{Float64}, 2)
	p[1] = zeros(n+1)
	p[2] = [a[:]; -1]
	D = Array(SparseMatrixCSC{Float64,Int}, 2)
	D[1] = sparse([covHalf zeros(n)]) 
	D[2] = sparse([3*RHalf zeros(n)]) 
	E = Array(SparseMatrixCSC{Float64,Int}, 2)
	w = Array(Vector{Float64}, 2)
	q = Array(Float64, 2)
	q[1] = -s
	q[2] = 0



	lx = zeros(n+1)  
	ux = ones(n+1)  
	if MIP
		f = zeros(n+1) 
		A = [speye(n+1); sparse([ones(1,n) 0]); spzeros(1,n+1)]
	 	B = [-speye(n+1); spzeros(1,n+1); sparse(ones(1,n+1))]		       
		b = [zeros(n+1); 1.0; k+1]
		sense = [['<' for i in [1:n+1]],'=','<'] 

		E[1] = zeros(n,n+1) 
		E[2] = zeros(n,n+1)
		w[1] = zeros(n+1)
		w[2] = zeros(n+1)

		lz = zeros(n+1)
		uz = ones(n+1)
	else
		f = Float64[]
		A = sparse([ones(1,n) 0])	
		B = spzeros(1,0)
		b = [1.0]
		sense = ['=']

		E[1] = zeros(n,0)
		E[2] = zeros(n,0)
		w[1] = zeros(0)
		w[2] = zeros(0)

		lz = Float64[]
		uz = Float64[]
	end
	misocp = MISOCPInput(c, f, A, B, b, sense, D, E, d, p, w, q, lx, ux, lz, uz)
	return misocp
end


function buildRobustMarkowitzPor(porfile::String,MIP=true)
	a, RHalf, covHalf=loadPorFileRobust(porfile)
	return buildRobustMarkowitzMISOCP(a, RHalf, covHalf,0.2,MIP)
end

function buildRobustMarkowitzOriginalPor(porfile::String,MIP=true)
	a, RHalf, covHalf=loadPorFileRobust(porfile)
	return buildRobustMarkowitzOriginalPorMISOCP(a, RHalf, covHalf,0.2,MIP)
end

function loadPorFileRobust(porfile::String)
	# Data = readdlm(porfile,' ')
	# n = int(Data[1,1])
	# covHalf= Data[3:(2+n),1:n]	
	# a = Data[2,1:n]
	# RHalf= Data[3+n:(2+2*n),1:n]

	# temp=diag(covHalf*covHalf')
	# @printf("Fraction of assets bellow risk threshold = %f\n", length(find(x->x<=0.2^2,temp))/n)
	# return a, RHalf, covHalf

	file = open(porfile,"r")
	n = int(readline(file))
	a = float(split(readline(file))[1:n])
	covHalf = zeros(n,n)
	for i in 1:n
		covHalf[i,:]=float(split(readline(file))[1:n])
	end
	RHalf = zeros(n,n)
	for i in 1:n
		RHalf[i,:]=float(split(readline(file))[1:n])
	end
	temp=diag(covHalf*covHalf')
	@printf("Fraction of assets bellow risk threshold = %f\n", length(find(x->x<=0.2^2,temp))/n)
	return a, RHalf, covHalf
end	




function buildModel(prob::MISOCPInput,whatSolver=MathProgBase.defaultQPsolver)
	model = Model(solver=whatSolver)

	nx = size(prob.A,2)
	nz = size(prob.B,2)
	m = size(prob.A,1)


	@defVar(model,prob.lx[i]<=x[i=1:nx]<=prob.ux[i])
	@defVar(model,prob.lz[i]<=z[i=1:nz]<=prob.uz[i],Int)
		
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
	for k in 1:nsoc
    	D = prob.D[k]' # transposed
		E = prob.E[k]' # transposed
        d = prob.d[k]
        p = prob.p[k]
		w = prob.w[k]
        q = prob.q[k]

		dim = length(d) # dimension of cone
        # y[k,1:dim] = D_k x +E_k z - d_k
        # y[k,0] = p_k^Tx +w_k^Tz - q_k
    	@defVar(model,y[k,i=1:dim])
    	@defVar(model,y0[k]>=0)
    	for i in 1:dim
            @addConstraint(model, sum{ D.nzval[idx]*x[D.rowval[idx]],
            	idx = D.colptr[i]:(D.colptr[i+1]-1)} + sum{ E.nzval[idx]*z[D.rowval[idx]],
            	idx = E.colptr[i]:(E.colptr[i+1]-1)}  - y[k,i] == d[i])
        end
        @addConstraint(model,  sum{ p[j]*x[j], j = 1:nx } + sum{ w[j]*z[j], j = 1:nz } - y0[k] == q)
        # macros don't yet accept quadratic terms
		addConstraint(model, sum([y[k,i]^2 for i=1:dim]) <= y0[k]^2)
	end

	return model, x, z
end







