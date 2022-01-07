import new_crank
import original_crank
import time
import matplotlib.pyplot as plt

def test_new_tridiag():
    time1 = time.time()
    for _ in range(200):
        new_crank.price_american_option_with_divs()
    time2 = time.time()
    print("New: %.4f seconds" % (time2-time1))

def test_old_method():
    time1 = time.time()
    for _ in range(200):
        original_crank.price_american_option_with_divs()
    time2 = time.time()
    print("Original: %.4f seconds" % (time2-time1))

# Tests 100 calculations of the original and new crank-nicolson method, measuring runtime
def test_crank():
    time1 = time.time()
    for _ in range(200):
        original_crank.price_american_option_with_divs()
    time2 = time.time()
    for _ in range(200):
        new_crank.price_american_option_with_divs()
    time3 = time.time()
    print("Original: %.4f seconds" % (time2-time1))
    print("New: %.4f seconds" % (time3-time2))

# Tests 2000 calculations each of the original and new implementations of the crank-nicolson method, measuring runtime every 100 calculations
def test_crank_plot():
    time_orig_start = time.time()
    original_times = [0]
    for i in range(20):
        for _ in range(100):
            original_crank.price_american_option_with_divs()
        original_times.append(time.time() - time_orig_start)
        print("loop %d complete, time= %.4f seconds" % (i+1, original_times[i+1]))
    time_new_start = time.time()
    new_times = [0]
    for i in range(20):
        for _ in range(100):
            new_crank.price_american_option_with_divs()
        new_times.append(time.time() - time_new_start)
        print("loop %d complete, time= %.4f seconds" % (i+1, new_times[i+1]))
    x = range(0, 2001, 100)
    plt.plot(x, original_times, "r^-", label="Original")
    plt.plot(x, new_times, "bo-", label="New")
    plt.legend()
    plt.title("Comparison of Crank-Nicolson Method Implementations\n(Ryzen 3700X, RTX 2080S, 32GB of RAM)")
    plt.xlabel("Number of calculations")
    plt.ylabel("Time (seconds)")
    plt.annotate("%.2f" % (original_times[-1]), xy=(x[-1], original_times[-1]), xytext=(x[-1], original_times[-1]))
    plt.annotate("%.2f" % (new_times[-1]), xy=(x[-1], new_times[-1]), xytext=(x[-1], new_times[-1]))
    plt.show() 
    plt.grid()
    plt.savefig("crank_plot.png")

if __name__ == "__main__":
    test_crank_plot()
    #test_new_tridiag()
