import numpy as np

# Fixed point algorithm
def fixedpt(x0,g,Nmax,Tol):
	
	N = 0
	xn1 = g(x0)
	xn = x0
	p = np.zeros((Nmax,1))
	# and np.abs(xn1-xn)>Tol
	while N<Nmax:
		p[N] = xn1
		xn = xn1
		xn1 = g(xn)
		N = N + 1
		print(N)
		return(xn1,p)
	
	
# Function to do fixed point on
def g(x): 
	return (-0.25*x+3)

# Initial guess and stopping parameters
x0 = 1
Nmax = 100
Tol = 1e-6

# Running fixed point
x,p = fixedpt(x0,g,Nmax,Tol)
print(p)
		
