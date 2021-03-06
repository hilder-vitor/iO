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

def mod_matrix(A, q):
    for i in xrange(A.nrows()):
        for j in xrange(A.ncols()):
            A[i, j] = my_modular(A[i, j], q)
    return A

    
# Return a 1xn matrix (row vector) with uniform random elements from Z_q
def random_vector_mod_q(n, q):
    return Matrix(1, n, [ZZ.random_element(-ceil(q/2), ceil(q/2), 'uniform') for i in xrange(n)])

def random_matrix_mod_q(nrows, ncols, q):
    A = matrix(nrows, ncols)
    for i in xrange(nrows):
        A[i] = [ZZ.random_element(0, q-1, 'uniform') for j in xrange(ncols)]
    return A
    
def gadget_matrix(n, q):
    len_g = ceil(log(q, 2))
    g = [2**i for i in xrange(0, len_g)]

    G = Matrix(n, n*len_g)

    for i in xrange(n):
        G[i] = [0] * len_g*i + g + [0] * len_g * (n-i-1)

    return G

def binary(u, n, q):
    l = ceil(log(q,2))
    resp = Matrix(n*l, 1)
    i = 0
    for u_i in u[:,0]:
        dig = u_i[0].digits(base=2, padto=l)
        for j in xrange(l):
            resp[i+j, 0] = dig[j]
        i += l
        
    return resp

def sample_SIS_matrix_and_trapdoors(n, m, q):
    l = ceil(log(q, 2))
    w = n*l
    if m <= w:
        return False
    m_bar = m - w

    A_bar = random_matrix_mod_q(n, m_bar, q)
    R = random_matrix_mod_q(m_bar, w, 2)
    G = gadget_matrix(n, q)

    A = block_matrix(1, 2, [A_bar, G - A_bar*R])

    return A, R # returning sampled matrix and its trapdoor

# Use trapdoor R to return small solution r to Ar = u (mod q) with A having format [A_bar | G - A_bar * R] 
def solve_ISIS_with_trapdoor(A, u, R, n, m, q):
    l = ceil(log(q, 2))
    w = n*l
    if m <= w:
        return False
    m_bar = m - w

    X = DiscreteGaussianDistributionIntegerSampler(sigma=ceil(sqrt(n*l)))
    rand_v = Matrix(m_bar, 1, [X() for i in xrange(m_bar)])
    #rand_v = random_matrix_mod_q(m_bar, 1, 10) # fix that bound: 10 is arbitrary
    
    r2 = binary((u - A[:,0:m_bar]*rand_v) % q, n, q)  # binary(u - A_bar*R)
    r1 = rand_v + R*r2

    return block_matrix(2, 1, [r1, r2])

# Use trapdoor R to return small solution r to Ar = U (mod q) with A having format [A_bar | G - A_bar * R] and U being a matrix
def solve_matrix_ISIS_with_trapdoor(A, U, R, n, m, q):
    solutions = []
    for j in xrange(U.ncols()):
        solutions.append(solve_ISIS_with_trapdoor(A, U[:, j], R, n, m, q))
    return block_matrix(1, U.ncols(), solutions)

def test_ISIS_solution(n, m, q):
    A, R = sample_SIS_matrix_and_trapdoors(n, m, q)
      
    U = random_matrix_mod_q(n, 3, q)
    
    for i in xrange(3):
        r = solve_matrix_ISIS_with_trapdoor(A, U, R, n, m, q)

        print "---- Is A*r = u % q ? ----"
        print A*r % q == U
        print r.nrows()
        print_mat(r)
        
    print "_________________________"

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

def run_tests():
    #####  LWE parameters #######
    n = 5
    q = next_prime(10*n)
    m = 10 * n * ceil(log(q, 2))
    alpha = 0.01  # not yet being used
    sigma = alpha * q
    #####  LWE parameters #######

    #X = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    print "-------- LWE (normal) instance---------"
    print "n =", n, " ....  m =", m, " ....  q =", q, " .... alpha =", alpha, " .... sigma =", sigma

    for i in xrange(3):
        test_ISIS_solution(n, m, q)
        
