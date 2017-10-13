from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler

k = 3
sec_param = 2

N = k+1 # XXX: check how to chose that value of N

def create_star_graph(k):
    """Create a graph to represent the k chains where the users will publish their individual
    values to get the shared secret.
    
    For each of the k users, there is a chain (one branch in the graph) with k vertices, each one connected to the next one in the chain
    and the last vertex of each chain is connected to the same vertex 0. Therefore, the returned graph is a "kind" of star graph with
    k*k + 1 vertices.

    For instance, for k=3, we have

    1,1 ->  1,2 -> 1,3 -------\
                               \
    2,1 ->  2,2 -> 2,3 --------> 0
                               /
    3,1 ->  3,2 -> 3,3 -------/
    
    """

    G = DiGraph(k*k + 1)
    
    for i in xrange(1,k+1):
        p = (i-1)*k+1
        for j in xrange(k-1):
            G.add_edge(p+j, p+j+1, str(i)+"-"+str(j))
            
    G.add_edge(p+k-1, 0, str(i)+"-0")
            
    return G


def generate_matrices_for_vertices(k, mmaps_params):
    As, trapdoors = instance_gen(params)
    C = {}
    for i in xrange(1,k+1):
        for l in xrange(N):
            t_i_l = sample(params)
            for i_prime in xrange(1,k+1):
                j = (i + i_prime) % k
                if 0 == j:   # edge from matrix A_i,k to A_0
                    source = i*k
                    dest = 0
                else:
                    source = (i-1)*k + j
                    dest = (i-1)*k + j + 1
                print "C[",i,", ",l,", ",i_prime,"]"
                C[i, l, i_prime] = encode(params, t_i_l, source, dest, As[source], As[dest], trapdoors[source])
    sources_each_chain = [As[(i-1)*k + 1] for i in xrange(k)]
    return C, sources_each_chain

def publish(public_encodings, distribution, source_user, dest_user):
    i_prime = source_user
    i = dest_user
    r_i = [distribution() for l in xrange(N)]
    D = public_encodings[i, 0, i_prime] * r_i[0]
    for l in xrange(1,N):
        D += public_encodings[i, l, i_prime] * r_i[l]
    return D




from graph_based_mmaps import params_gen
from graph_based_mmaps import instance_gen
from graph_based_mmaps import sample
from graph_based_mmaps import encode
from graph_based_mmaps import zero_test
from graph_based_mmaps import extract


G = create_star_graph(k)
params = params_gen(sec_param, G)

C, sources_each_chain = generate_matrices_for_vertices(k, params)

print C
print sources_each_chain

X=params['Chi']

D_0_1 = publish(C, X, 1, 0)
print D_0_1
D_1_1 = publish(C, X, 1, 1)
print D_1_1
D_2_1 = publish(C, X, 1, 2)
print D_2_1
D_3_1 = publish(C, X, 1, 3)
print D_3_1

D_0_2 = publish(C, X, 2, 0)
print D_0_2
D_1_2 = publish(C, X, 2, 1)
print D_1_2
D_2_2 = publish(C, X, 2, 2)
print D_2_2
D_3_2 = publish(C, X, 2, 3)
print D_3_2

D_0_3 = publish(C, X, 3, 0)
print D_0_3
D_1_3 = publish(C, X, 3, 1)
print D_1_3
D_2_3 = publish(C, X, 3, 2)
print D_2_3
D_3_3 = publish(C, X, 3, 3)
print D_3_3


A1 = sources_each_chain[0]
combined = D_2_1 * D_0_1 * D_1_1 * A1

secret = extract(params, combined, A1)

print secret
