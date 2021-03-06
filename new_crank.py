from numpy import *
from math import *
import matplotlib.pyplot as plt
import time

from numpy.linalg import solve
from numba import jit, cuda, f8, njit
import warnings
warnings.filterwarnings("ignore")

#JAY -- using Numba, the JIT compiler stores the function in a compiled binary that dramatically speeds up runtime
@jit(f8[:] (f8[:],f8[:],f8[:],f8[:] ), cache=True, nopython=True)
def solve_tridiagonal(top, mid, bottom, b):
    '''
    JAY -- a fast function for solving tridiagonal matrices, based on Thomas' Algorithm. 
    '''
    height = len(b)

    for i in range(1, height):
        m = top[i-1]/mid[i-1]
        mid[i] = mid[i] - m*bottom[i-1] 
        b[i] = b[i] - m*b[i-1]
        	    
    x = mid
    x[-1] = b[-1]/mid[-1]

    for i in range(height-2, -1, -1):
        x[i] = (b[i]-bottom[i]*x[i+1])/mid[i]

    return x

@jit(cache=True)
def price_american_option_with_divs():

    #JAY -- Initialize variables -- can tweak these to compare against known pricings
    S = 4
    Smin = 0.4
    Smax = 150
    E = 8
    r = 0.1
    d = 0.08
    sigma = 0.4
    T = 1
    t = 0
    X = log(S/E)
    k = r/(0.5*sigma**2)
    kd = (r-d)/(0.5*sigma**2)

    M = 100
    Nminus = -100
    Nplus = 100
    N = Nplus-Nminus
    u = zeros((N+1, M+1))
    dt = (0.5*sigma**2*T)/M
    Xzero = log(Smin/E)
    Xmax = log(Smax/E)
    dX = (Xmax-Xzero)/N
    d_plot_X = (Smax - Smin)/N
    alpha = dt/(dX*dX)

    #JAY -- Used for measuring how long the program will take
    start_time = time.time()

    Xmesh = [Xzero + dX*i for i in range(N+1)]
    X_plot_mesh = [Smin + d_plot_X*i for i in range(N+1)]

    g = zeros((N+1, M+1))

    #JAY -- Initialize the payoff matrix
    for n in range(1, N+2):
        for m in range(2, M+2):
            g[n-1, m-1] = exp((0.25*(kd-1)**2+k)*((m-1)*dt))*(max((exp(0.5 * (kd+1)*(n+Nminus-1)*dX)-exp(0.5*(kd-1)*(n+Nminus-1)*dX)), 0))
        g[n-1, 0] = max((exp(0.5*(kd+1)*(n+Nminus-1)*dX) -
                      exp(0.5*(kd-1)*(n+Nminus-1)*dX)), 0)
    g[0] = [0]*(M+1)

    #JAY -- Boundary conditions: note that the right side is currently unbounded (Potential Optimization?)
    u[N] = g[N]
    u[:, 0] = g[:, 0]

    a = -alpha
    b = 1+2*alpha
    c = -alpha

    zmatrix = zeros((N+1, M+1))

    #JAY -- the tridiagonal matrix (as shown in the original paper), split into each diagonal (note that the top and bottom diagonal are identical)
    top_and_bottom = [-0.5*alpha] * (N-2)
    mid = [1+alpha] * (N-1)


    #JAY -- The main part of the logic; iterates starting at t=0 until the expiry time, determining the value of the next time step using the crank-nicolson method of finite differences
    for p in range(M):
        temp = zeros(N-1)
        temp[0] = a*g[0, p+1]
        temp[-1] = c*g[N, p+1]
        zmatrix[1:N, p] = (1-alpha)*u[1:N, p]+0.5 * \
            alpha*(u[2:N+1, p]+u[:N-1, p])
        RHS = zmatrix[1:N, p] - temp
        b = RHS

        #JAY -- <-MAIN OPTIMISATION-> The original paper used an iterative SOR algorithm which was inefficient, numpy's matrix solver for square matrices allows for parallel processing and is therefore much faster, and also does not require multiple iterations. This can be further optimised with a tridiagonal solver, which can solve the system of equations in linear time
        x = solve_tridiagonal(array(top_and_bottom), array(mid), array(top_and_bottom), array(b))
        u[1:N, p+1] = x

    #JAY -- This is the final value of the option, which is the interpolated value of the last time step
    uresult = interp(X, Xmesh, u[:, M])

    #JAY -- End result is obtained using the Black-Scholes formula
    return (E*E**(0.5*(kd-1))*S**(-0.5*(kd-1))*exp((-(1/4)*(kd-1)**2-k)*0.5*sigma**2*T)*uresult, time.time()-start_time)


if __name__ == "__main__":
    print("Price of American Option: %.4f\nTime Taken to calculate: %.4f seconds" %
          price_american_option_with_divs())
