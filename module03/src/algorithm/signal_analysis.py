

def analyse_homo_signals(data):
    # compute homology signals of interaction sites based on site adjacency and peptide overlap
    #
    # input data: pd.DataFrame
    # return data: pd.DataFrame

    data["homo_adjacency"] = data.apply(lambda x: compute_interaction_adj(x), axis=1)
    data["homo_int_overl"] = data.apply(lambda x: compute_interaction_overlap(x), axis=1)
    data["homo_pep_overl"] = data.homo_int_overl > 0
    return data


def compute_interaction_adj(data_row):
    # compute interactions site adjacency, represented by value between 0 (residues far away on sequence) and 1
    # (residues are the same)
    #
    # input data_row: pd.Series
    # return compute_interaction_dist: float

    return 1 - (abs(int(data_row["pos_a"]) - int(data_row["pos_b"])) / len(data_row["seq_a"]))


def compute_interaction_overlap(data_row):
    # compute peptide overlap between/including interacting residues, represented by value between 0
    # (no peptide overlap between/including interacting residues) and 1 (both interacting residues are in both peptides)
    #
    # input data_row: pd.Series
    # return compute_interaction_overlap: float

    if data_row["pos_a"] == data_row["pos_b"]:
        return 1.0
    else:
        site1 = 'a' if data_row["pos_a"] < data_row["pos_b"] else 'b'
        site2 = 'a' if data_row["pos_a"] > data_row["pos_b"] else 'b'
        seq, pep_a, pep_b, pos_a, pos_b = (data_row[f"seq_{site1}"], data_row[f"pep_{site1}"], data_row[f"pep_{site2}"],
                                           int(data_row[f"pos_{site1}"]) - 1, int(data_row[f"pos_{site2}"]) - 1)

        # save indices of residues in peptides between/including interacting residues
        seq_a_inds = [ind for ind in range(seq.find(pep_a), seq.find(pep_a) + len(pep_a)) if ind >= pos_a]
        seq_b_inds = [ind for ind in range(seq.find(pep_b), seq.find(pep_b) + len(pep_b)) if ind <= pos_b]

        # compute intersect and union of index lists
        seq_intersect = [ind for ind in seq_a_inds if ind in seq_b_inds]
        seq_union = seq_a_inds.copy()
        seq_union.extend([ind for ind in seq_b_inds if ind not in seq_a_inds])

        if len(seq_intersect) == 0:
            return 0.0
        else:
            return len(seq_intersect) / len(seq_union)