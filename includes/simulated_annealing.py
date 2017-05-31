# The various temparature functions for simulated annealing

# Temperature stays maximal, which means that simulated annealing with oneT is
# basically the same as our greedy approach
def gen_oneT():
    def oneT(i):
        return 1
    return oneT

# We trya linear cooldown of the Temperature
def gen_linearT(maxi):
    def linearT(i):
        return i/maxi
    return linearT

# We try a coolingdown that starts small, but then gets cold very quickly
def gen_exponentialT(maxi, b):
    g = (1/b)**(1/maxi)
    def exponentialT(i):
        return b*g**i
    return exponentialT

# f:         y0 + (y1 - y0)
#    -----------------------------
#                        x0+x1
#     (1 + exp(-v *(i -  ----- )
#                          2
#
# We try different approaches to calculate v, the
# average steepness of the slope

# We calculated v algebraicly
def gen_sigmoidT_mathv(maxi, eps=10e-4):
    x1 = maxi/2
    v = 2*np.log((1/eps) - 1) / maxi
    def sigmoidT(i):
        return 1 / (1 + np.exp(-v * (i - x1)))
    return sigmoidT

# We approach v by checking if it meets certain requirements
def gen_sigmoidT_tryv(maxi, eps=10e-6):
    x1 = maxi/2
    v = 1
    def sigmoidT(i):
        return 1 / (1 + np.exp(-v * (i - x1)))

    while not np.abs(    sigmoidT(0   )) < eps or\
          not np.abs(1 - sigmoidT(maxi)) < eps:
        v /= 2
    return sigmoidT
