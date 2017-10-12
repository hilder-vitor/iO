
def matrix_norm(A):
    max_val = 0
    for i in xrange(A.nrows()):
        for j in xrange(A.ncols()):
            if abs(A[i,j]) > max_val:
                max_val = abs(A[i,j])            

    return max_val

def create_digraph(number_of_vertices, list_of_edges):
    G = DiGraph(number_of_vertices)
    for e in list_of_edges:
        G.add_edge(e)
    return G

def diameter(G):
    '''The sage's method Graph.diameter returns infinity if the digraph is not strongly connected.
    This function returns -1 if the graph has no edge or the maximum distance between all connected vertices.'''
    dists = G.shortest_path_all_pairs()[0] # distances ignoring +Infinity (non connected vertices)
    diameter = -1
    for h in dists.values():
        d = max(h.values())
        if d > diameter:
            diameter = d
    return diameter

def sample_matrix(nrow, ncol, dist):
    return Matrix(nrow, ncol, [[dist() for j in xrange(ncol)] for i in xrange(nrow)])


def params_gen(sec_param, graph):
    d = diameter(graph)
    n = ceil(d * sec_param * log(d * sec_param, 2))
    q = next_prime((d * sec_param) ^ d)
    m = n * d * ceil(log(q, 2))
    s = sqrt(n)
    sigma = sqrt(n * (d + 1) * ceil(log(q, 2))) / 4
    t = floor(log(q, 2) / 4) - 1
    X = DiscreteGaussianDistributionIntegerSampler(sigma=sigma)
    return {'lambda':sec_param, 'd':d, 'n':n, 'q':q, 'm':m, 's':s, 'sigma':sigma, 't':t, 'G':graph, 'Chi':X}

from lattice_trapdoor_SIS import sample_SIS_matrix_and_trapdoors
from lattice_trapdoor_SIS import solve_matrix_ISIS_with_trapdoor
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
def instance_gen(params):
    matrices = []
    trapdoors = []
    for i in xrange(params['G'].order()):
        A, R = sample_SIS_matrix_and_trapdoors(params['n'], params['m'], params['q'])
        matrices.append(A)
        trapdoors.append(R)

    return matrices, trapdoors

def sample(params):
    return sample_matrix(params['n'], params['n'], params['Chi'])

def encode(params, S, u, v, Au, Av, trapdoor_u):
    E = sample_matrix(params['n'], params['m'], params['Chi'])
    D = solve_matrix_ISIS_with_trapdoor(Au, S*Av + E, trapdoor_u, params['n'], params['m'], params['q'])
    return D

# Test if D encodes zero relative to some path begining in u
from lattice_trapdoor_SIS import mod_matrix
def zero_test(params, D, Au):
    if matrix_norm(mod_matrix(Au*D, params['q'])) < params['q'] / (2^(params['t']+1)):
        return True
    return False


import hashlib
def extract(params, D, Au):
    max_value = matrix_norm((Au*D) % params['q'])

    extracted_bits = str(max_value // (2^(params['t']+1)))

    hex_extracted_val = hashlib.sha1(extracted_bits).hexdigest()

    return ZZ('0x' + hex_extracted_val)
   


########################################################################
####        MAIN
def test_graph_based_mmap():
    security_parameter = 9
    G = create_digraph(4, [[1,2], [1,3], [2,0]])

    params = params_gen(security_parameter, G)
    print params

    matrices, trapdoors = instance_gen(params)

    S = sample(params)

    print "D1 = encode(params, S, 1, 3, A_1, A_3, trapdoors[1])"
    D1 = encode(params, S, 1, 3, matrices[1], matrices[3], trapdoors[1])

    print "extract(params, D1, matrices[1])"
    extracted_1 = extract(params, D1, matrices[1])
    print extracted_1

    print "D2 = encode(params, S, 1, 3, A_1, A_3, trapdoors[1])"
    D2 = encode(params, S, 1, 3, matrices[1], matrices[3], trapdoors[1])

    print "extract(params, D2, matrices[1])"
    extracted_2 = extract(params, D2, matrices[1])
    print extracted_2

    print "Do the two random extracted values are equal?"
    print extracted_1 == extracted_2

    D3 = D1 - D2
    print "Does D1 - D2 encodes zero?"
    print zero_test(params, D3, matrices[1])
