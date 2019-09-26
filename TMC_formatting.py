import dendropy


def map_cell_ids_for_sagi(rtd):
    """
    Translate cell names to 0..n TMC input convention
    """
    cell_id_map_for_sagi = dict()
    for i, cell_id in enumerate(sorted(rtd.keys())):
        cell_id_map_for_sagi[cell_id] = i
    return cell_id_map_for_sagi


def format_triplet(triplet, pair, score, cell_id_map_for_sagi, print_scores=True):
    """
    <index of leaf1>, <index of leaf2> | <index of leaf3>[:<weight>] when leaf1 and leaf2 are brothers and leaf3 is a cousin?
    """

    if pair is None:
        return ''
    brother_a, brother_b = pair
    triplets_labels = triplet
    cousins = triplets_labels - {brother_a, brother_b}
    cousin, = cousins  # implicit assertion that len(cousins) == 1
    triplet_string = '{},{}|{}'.format(cell_id_map_for_sagi[brother_a], cell_id_map_for_sagi[brother_b], cell_id_map_for_sagi[cousin])
    assert len(triplet_string) <= 15
    scores_string=''
    if print_scores:
        scores_string=':{:.4f}'.format(score)
    delimiter_string=' '
    return triplet_string + scores_string + delimiter_string


def convert_tree_with_cell_id_map(taxon_name_space, cell_id_map_for_sagi):
    assert len(cell_id_map_for_sagi.keys()) == len(set(cell_id_map_for_sagi.values()))
    reverse_cell_id_map_for_sagi = {v: k for k, v in cell_id_map_for_sagi.items()}
    label_taxon_map = taxon_name_space.label_taxon_map()
    for l in label_taxon_map:
        original_label = reverse_cell_id_map_for_sagi[int(l)]
        taxa = label_taxon_map[l]
        taxa.label = original_label


def convert_names_in_sagis_newick(indexes_as_labels_tree, labeled_tree, cell_id_map_for_sagi):
    """
    Convert 0..n labels in the TMC resulting newick back to the original labels
    """
    tns_sagi = dendropy.TaxonNamespace()
    indexes_tree = dendropy.Tree.get_from_path(
        indexes_as_labels_tree,
        "newick",
        taxon_namespace=tns_sagi)
    convert_tree_with_cell_id_map(tns_sagi, cell_id_map_for_sagi)
    indexes_tree.write_to_path(labeled_tree, 'newick', suppress_rooting=True)