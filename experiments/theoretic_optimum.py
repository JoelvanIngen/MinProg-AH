_bond_values = {
    frozenset({'H', 'H'}): -1,
    frozenset({'C', 'C'}): -5,
    frozenset({'C', 'H'}): -1,
}


def theoretic_optimum(protein_str: str):
    """
    Calculates the theoretic optimum for a given protein string by iterating over the string and taking
     the two highest possible bonds a node can make.
     Author: Jarec

    :param protein_str: the protein string
    :return: score: the theoretic optimum score
    """
    score = 0
    for i, letter in enumerate(protein_str[:-3], start=3):
        possible_bonds = [_bond_values.get(frozenset({letter, l}), 0) for l in protein_str[i::2]]
        possible_bonds.sort()
        score += sum(possible_bonds[:2])
    return score


if __name__ == '__main__':
    protein = "HHCPCHPCHPCHHHCPCH"
    print(theoretic_optimum(protein))
