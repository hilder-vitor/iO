def my_modular(x, q):
    if 0 <= x and x <= q-1:
        if x > q/2:
            return x - q
        else:
            return x
    else:
        new_x = x % q # usual modular operation
        return my_modular(new_x, q)

def mod_vec(v, q):
    return Matrix(1, v.ncols(), [my_modular(v[0,i], q) for i in xrange(v.ncols())])
    
    
def gadget_matrix(n, q):
    len_g = ceil(log(q, 2))
    g = [2**i for i in xrange(0, len_g)]

    G = Matrix(n, n*len_g)

    for i in xrange(n):
        G[i] = [0] * len_g*i + g + [0] * len_g * (n-i-1)

    return G


def trap_door_for_gadget_matrix(n, q):
    bin_q = q.digits(2)
    log_q = ceil(log(q, 2))

    Tg = Matrix(log_q, log_q)

    for i in xrange(1,log_q-1):
        Tg[i,i-1] = -1
        Tg[i,i] = 2
        Tg[i,log_q-1] = bin_q[i]

    Tg[0,0] = 2
    Tg[0,log_q-1] = bin_q[0]

    Tg[log_q-1,log_q-2] = -1
    Tg[log_q-1,log_q-1] = bin_q[log_q-1]
    
    In = matrix.identity(n)
    TG = In.tensor_product(Tg)

    return TG

def encode_LWE_given_matrix(G, s, error_dist, n, m, q):
    X = error_dist
    e = Matrix(1, m, [X() for i in xrange(m)])
    return (G, ((s*G + e) % q))


def solve_LWE_gadget_with_trapdoor(G, b, T, n, m, q):
    invT = T.inverse()
    sG = mod_vec(b - mod_vec(b*T, q) * invT, q)
    return Matrix(1, n, [sG[0, i] for i in xrange(0, m, m/n)])



from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

#####  LWE parameters #######
n = 3
q = next_prime(32)
m = n * ceil(log(q, 2))
alpha = 0.02
sigma = alpha * q
#####  LWE parameters #######

X = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
G = gadget_matrix(n, q)
T = trap_door_for_gadget_matrix(n, q)

s = [ZZ.random_element(-ceil(q/2), ceil(q/2), 'uniform') for i in xrange(n)]
s = Matrix(1, n, s)

G, b = encode_LWE_given_matrix(G, s, X, n, n*(ceil(log(q,2))), q)

s_found = solve_LWE_gadget_with_trapdoor(G, b, T, n, m, q)


print "-------- LWE instance---------"
print "n =", n, " ....  m =", m, " ....  q =", q, " .... alpha =", alpha, " .... sigma =", sigma
print "G ="
print G
print "b =", b
print "s =", s

print "------- Trapdoor ------------"
print T

print "------- s found ------------"
print s_found

print "------- Is this solution OK? ------------"
print s == s_found
