def print_vec(v):
    for x in v:
        print x,
    print ""

def print_mat(A):
    for v in A:
        print_vec(v)


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
    
# Return a 1xn matrix (row vector) with uniform random elements from Z_q
def random_vec_mod_q(n, q):
    return Matrix(1, n, [ZZ.random_element(-ceil(q/2), ceil(q/2), 'uniform') for i in xrange(n)])

def random_matrix_mod_q(nrows, ncols, q):
    A = matrix(nrows, ncols)
    for i in xrange(nrows):
        A[i] = [ZZ.random_element(-ceil(q/2), ceil(q/2), 'uniform') for j in xrange(ncols)]
    return A
    
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

def encode_LWE_given_matrix(A, s, error_dist, n, m, q):
    X = error_dist
    e = Matrix(1, m, [X() for i in xrange(m)])
    # FIXME: this e is returned for debugging porpuse
    return (A, mod_vec(s*A + e, q), e)


def sample_LWE_matrix_and_trapdoors(n, m, q):
    m0 = m + 2*n
    A0 = random_matrix_mod_q(n, m0, q)
    R0 = random_matrix_mod_q(m0, m, 2)

    X = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    G = gadget_matrix(n, q)
    print "---G:"
    print G
    T = trap_door_for_gadget_matrix(n, q)
    print "---Original trapdoor:"
    print T

    A = block_matrix(1, 2, [A0, G - A0*R0])
    print A.ncols()
    R = block_matrix(2, 1, [[R0], [matrix.identity(m, m)]])

    return A,R,T

def solve_LWE_gadget_with_trapdoor(G, b, T, n, m, q):
    invT = T.inverse()
    sG = mod_vec(b - mod_vec(b*T, q) * invT, q)
    return Matrix(1, n, [sG[0, i] for i in xrange(0, m, m/n)])

def solve_LWE_with_trapdoors(A, b, R, T, n, m, q):
    # XXX: This will not work if b*R*T  has some entry bigger than q/2
    return solve_LWE_gadget_with_trapdoor(A*R, mod_vec(b*R,q), T, n, m, q)



from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

#####  LWE parameters #######
n = 5
q = next_prime(10*n*n)
m = n * ceil(log(q, 2))
alpha = 0.01  # XXX: when alpha is bigger, sigma becomes "big" and b*R*T mod q isn't equal to b*R*T
sigma = alpha * q
#####  LWE parameters #######

X = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
G = gadget_matrix(n, q)
T = trap_door_for_gadget_matrix(n, q)

s = random_vec_mod_q(n, q)

G, b, E = encode_LWE_given_matrix(G, s, X, n, m, q)

s_found = solve_LWE_gadget_with_trapdoor(G, b, T, n, m, q)

A,R,T = sample_LWE_matrix_and_trapdoors(n, m, q)
print A
print R
print T

s = random_vec_mod_q(n, q)
# FIXME: this E is returned for debugging porpuse
A, b, E = encode_LWE_given_matrix(A, s, X, n, 2*m+2*n, q)

print "E ="
print_vec(E)
print "E*R = "
print_vec(E*R)
s_found = solve_LWE_with_trapdoors(A, b, R, T, n, m, q)
print "E*R*T = "
print_vec(E*R*T)

print "-------- LWE (normal) instance---------"
print "n =", n, " ....  m =", m, " ....  q =", q, " .... alpha =", alpha, " .... sigma =", sigma
print "A ="
print A
print "b =", b
print "s =", s

print "------- Trapdoor ------------"
print T

print "------- s found ------------"
print s_found

print "------- Is this solution OK? ------------"
print s == s_found



