import numpy as np

# to calculate a European call, just use American call, they should be equal
def op():
    EA = 1
    # CP = int(input("Enter 0 for call, 1 for put: "))
    # S = float(input("Current price: "))
    # K = float(input("Strike: "))
    # T = float(input("Expiration: "))
    # sig = float(input("Volatility: "))
    # r = float(input("Risk-free interest rate: "))
    # q = float(input("Continuous dividend rate: "))
    # M = int(input("Number of price steps: "))
    # N = int(input("Number of time steps: "))
    CP = 0
    S = 15
    K = 8
    T = 1
    sig = 0.4
    r = 0.1
    q = 0.08
    M = 200
    N = 200

    # price and time step widths with mult*S as the max stock price
    mult = 5

    # change M to make it easier to recover option price at the end
    M -= M % mult
    Smax = mult*S

    s = Smax/M
    t = T/N

    # initialize the grid
    g = np.zeros((M+1, N+1))

    for i in range(N+1):
        g[0][i] += Smax-K*np.exp(-(r-q)*(T-t*i))

    for i in range(1, M):
        g[i][N] = max((M-i)*s-K, 0)

    g[M - M//mult][0] = S
    g[M - M//mult][N] = max(S-K, 0)

    # define the coefficients to be used in solving difference equations
    # note that we need to switch the j indices around: in my head I'm thinking
    # of lower indices corresponding to higher up on the y-axis, just like
    # in the x-y plane
    def a(j):
        return -0.5*(r-q)*(M-j)*t - 0.5*sig**2*(M-j)**2*t

    def b(j):
        return 1 + sig**2*(M-j)**2*t + r*t

    def c(j):
        return 0.5*(r-q)*(M-j)*t - 0.5*sig**2*(M-j)**2*t

    # work backwards from the righthand side, updating grid column-by-column
    for i in reversed(range(1,N+1)):
        A = np.identity(M+1)
        d = g[:,i]

        for j in range(1, M):
            A[j][j-1] = a(j)
            A[j][j] = b(j)
            A[j][j+1] = c(j)

        f = np.linalg.solve(A, d)
        f = [ ( (1-CP)*max(f[k], (M-k)*s - K) + CP*(max(f[k], EA*(K-(M-k)*s))) ) for k in range(M+1)]

        # print out A, d, and f used in Af = d to use in debugging
        # print("\n\n", A, "\n\n", d, "\n\n", f)

        g[1:M,i-1] = f[1:M]

    # np.set_printoptions(precision=3)

    print(g[M-M//mult][0])

if __name__ == "__main__":
    op()