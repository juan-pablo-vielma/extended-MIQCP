function covCor(x,shrink=-1)

# x (t*n): t iid observations on n random variables
# sigma (n*n): invertible covariance matrix estimator
#
# Shrinks towards constant correlation matrix
# if shrink is specified, then this constant is used for shrinkage

# The notation follows Ledoit and Wolf (2003, 2004)
# This version 02/2010

# de-mean returns

(t,n)=size(x)
meanx=mean(x,1)
for i in 1:t
	x[i,:] -= meanx;
end

# compute sample covariance matrix
sample=(1/t).*(x'*x)

# compute prior
var=diag(sample)
sqrtvar=sqrt(var)

rBar = 0

prior = zeros(n,n)
for i in 1:n
	for j in 1:n
		prior[i,j] = sqrtvar[i]*sqrtvar[j]
		rBar += ( (sample[i,j]) / prior[i,j] );
	end
end

rBar = ( rBar - n ) / (n*(n-1))

prior *= rBar

for i = 1:n
	prior[i,i] = var[i];
end

if (shrink == -1) # compute shrinkage parameters and constant
			       
  # what we call pi-hat
  y=x.^2
  phiMat=y'*y/t - 2*(x'*x).*sample/t + sample.^2
  phi=sum(phiMat)
  
  # what we call rho-hat
  term1=((x.^3)'*x)/t
  thetaMat = zeros(n,n)
  help = x'*x/t
  helpDiag=diag(help)
  for i = 1:n
  	for j = 1:n
		thetaMat[i,j] = term1[i,j]
		thetaMat[i,j] -= helpDiag[i]*sample[i,j]
		thetaMat[i,j] -= help[i,j]*var[i]
		thetaMat[i,j] += var[i]*sample[i,j]
	end
  end
  for i = 1:n
  	thetaMat[i,i] = 0
  end
  rho=sum(diag(phiMat))+rBar*sum(((1./sqrtvar)*sqrtvar').*thetaMat);
  
  # what we call gamma-hat
  gamma=vecnorm(sample-prior)^2; #normfro(sample-prior)^2;
  
  # compute shrinkage constant
  kappa=(phi-rho)/gamma;
  shrinkage=max(0,min(1,kappa/t));
  
else # use specified constant
  shrinkage=shrink;
end

# compute the estimator
sigma=shrinkage*prior+(1-shrinkage)*sample;

return sigma, shrinkage

end

#x = float([ 0 0 1 0 1; 1 1 1 0 1; 0 1 0 1 0 ])
#(sigma, shrinkage) = covCor(x,-1)
#println(sigma)
