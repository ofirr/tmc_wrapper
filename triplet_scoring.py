import math
from itertools import combinations
from .utils import subdict

class NoIntersectionLociException(Exception):
    pass


def uri_prep(root, cella, cellb):
    """
    Calculate presets for triplet scoring
    """
    loci = cella.keys() & cellb.keys() & root.keys()  # safer, assumes imperfect root
    if not loci:
        raise NoIntersectionLociException
    n = float(len(loci))
    cella = subdict(cella, loci)
    cellb = subdict(cellb, loci)
    root = subdict(root, loci)
    sab = float(len(cella.items() & cellb.items() - root.items()))
    sa = float(len(cella.items() - root.items()))
    sb = float(len(cellb.items() - root.items()))
    return n, sab, sa, sb


def uri10(root, cella, cellb):
    """
    Approximate the distance between the pair MRCA and the root
    """
    n, sab, sa, sb = uri_prep(root, cella, cellb)
    if sa and sb:
        return (sab - (sa+sb)/n)/math.sqrt(sa*sb/n)
    raise NoIntersectionLociException


def choose_best_pair(trip, d):
    pair_scores = []
    for pair in combinations(trip, 2):
        cell_a, cell_b = pair
        pair_score = uri10(d['root'], d[cell_a], d[cell_b])
        pair_scores.append((pair, pair_score))
    pair = sorted(pair_scores, key=lambda x:x[1], reverse=True)  # sort by pair score (higher is better)
    return pair[0][0], pair[0][1]-pair[1][1]