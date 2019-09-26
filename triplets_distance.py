import itertools
from dendropy import SeedNodeDeletionException
from .utils import memory_expensive_random_choose

def check_triplet(tree, nodes_triplet, ref_close_pair):
    try:
        rec_triplet = tree.extract_tree_with_taxa(nodes_triplet)
    except SeedNodeDeletionException:
        return None
    if len(rec_triplet.leaf_nodes()) < 3:
        return None
    rec_pdm = tree.phylogenetic_distance_matrix()
    pairs_in_triplet = list(itertools.combinations(nodes_triplet, 2))
    sorted_pairs_and_edge_counts = [(pair, rec_pdm.path_edge_count(*pair)) for pair in pairs_in_triplet]
    sorted_pairs_by_edge_counts = sorted(sorted_pairs_and_edge_counts, key=lambda x: x[1])
    if sorted_pairs_by_edge_counts[0][1] == sorted_pairs_by_edge_counts[1][1]:  # tie, this can happen with non-binary furcations
        return None
    rec_close_pair = set(sorted_pairs_by_edge_counts[0][0])  # you get the closest
    return ref_close_pair == rec_close_pair



def challenge_triplets_generator(ref_tree, ndm, pdm, n=1000, excluded_nodes=None):
    nodes = {n.taxon for n in ref_tree.leaf_nodes()}
    if excluded_nodes is not None:
        nodes = nodes - set(excluded_nodes)

    for nodes_triplet in memory_expensive_random_choose(nodes, 3, n=n):
        pairs_in_triplet = list(itertools.combinations(nodes_triplet, 2))
        ref_close_pair = set(sorted(pairs_in_triplet, key=lambda x: pdm.path_edge_count(x[0],x[1]))[0])
        yield nodes_triplet, ref_close_pair


def triplets_score(tree_in_question, reference_tree, n=1000, excluded_nodes=None):
    results_by_difficulty = {True: 0, False: 0, None: 0}
    ndm = reference_tree.node_distance_matrix()
    pdm = reference_tree.phylogenetic_distance_matrix()
    for nodes_triplet, ref_close_pair in list(
        challenge_triplets_generator(
            reference_tree, 
            ndm, 
            pdm,
            n=n, 
            excluded_nodes=excluded_nodes)):
        check_status = check_triplet(tree_in_question, nodes_triplet, ref_close_pair)
        results_by_difficulty[check_status] += 1
    return results_by_difficulty