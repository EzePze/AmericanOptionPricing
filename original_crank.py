from numpy import *
from math import *
import matplotlib.pyplot as plt
import time


def price_american_option_with_divs():
    #S = int(input("Current price: "))
    S = 15
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

    tol = 0.001
    omega = 1.2

    # fig, ax = plt.subplots()
    # ax.grid()
    # plot_step = 25

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

    start_time = time.time()

    #Xmesh = [Xzero, dX, Xmax]
    #Xmesh = range(Xzero, Xmax, dX)
    Xmesh = [Xzero + dX*i for i in range(N+1)]
    X_plot_mesh = [Smin + d_plot_X*i for i in range(N+1)]

    g = zeros((N+1, M+1))

    for n in range(N+1):
        for m in range(1, M+1):
            g[n, m] = exp((0.25*(kd-1)**2+k)*((m-1)*dt))*(max((exp(0.5 *
                                                                   (kd+1)*(n+Nminus-1)*dX)-exp(0.5*(kd-1)*(n+Nminus-1)*dX)), 0))
        g[n, 0] = max((exp(0.5*(kd+1)*(n+Nminus-1)*dX) -
                      exp(0.5*(kd-1)*(n+Nminus-1)*dX)), 0)
    g[0] = [0]*(M+1)

    u[N] = g[N]
    u[:, 0] = g[:, 0]

    a = -alpha
    b = 1+2*alpha
    c = -alpha

    zmatrix = zeros((N+1, M+1))

    cmatrix = diagflat([-0.5*alpha for _ in range(N-2)], -1) +\
        diagflat([1+alpha for _ in range(N-1)]) +\
        diagflat([-0.5*alpha for _ in range(N-2)], 1)

    for p in range(M):
        temp = zeros(N-1)
        temp[0] = a*g[0, p+1]
        temp[-1] = c*g[N, p+1]
        zmatrix[1:N, p] = (1-alpha)*u[1:N, p]+0.5 * \
            alpha*(u[2:N+1, p]+u[:N-1, p])
        RHS = zmatrix[1:N, p] - temp
        b = RHS
        #print(u[1:N, p].shape)
        #x = max(u[1:N][p], g[1:N][p])
        x = array([i if i > j else j for i, j in zip(u[1:N, p], g[1:N, p+1])])
        # print(x.shape)
        xold = 1000*x
        n = len(x)
        while linalg.norm(xold - x) > tol:
            xold = x
            for i in range(n):
                if i == 0:
                    z = (b[i]+0.5*alpha*x[i+1])/(1+alpha)
                    x[i] = max(omega*z + (1-omega)*xold[i], g[i][p])

                elif i == n-1:
                    z = (b[i]+0.5*alpha*x[i-1])/(1+alpha)
                    x[i] = max(omega*z + (1-omega)*xold[i], g[i][p])

                else:
                    z = (b[i]+0.5*alpha*(x[i-1]+x[i+1]))/(1+alpha)
                    x[i] = max(omega*z + (1-omega)*xold[i], g[i][p])

        x = linalg.solve(cmatrix, b)
        u[1:N, p+1] = x
        # if p % plot_step == 0:
        # ax.plot(X_plot_mesh, u[:, M], label="t = %.2f" % (p*dt))
        # ax.set(xlabel='S', ylabel='u(X,t)', title='u(X,t)')
        # fig.savefig('crank-nicholson-calc.png')
        # plt.show()

    uresult = interp(X, Xmesh, u[:, M])

    return (E*E**(0.5*(kd-1))*S**(-0.5*(kd-1))*exp((-(1/4)*(kd-1)**2-k)*0.5*sigma**2*T)*uresult, time.time()-start_time)


if __name__ == "__main__":
    print("Price of American Option: %.4f\nTime Taken to calculate: %.4f seconds" %
          price_american_option_with_divs())
