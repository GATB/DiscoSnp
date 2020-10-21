import sys


def gfa_to_fa(gfa_file_name: str):
    """Prints genomic sequences from GFA nodes

    Args:
        gfa_file_name (str): File name of the input GFA
    """

    with open(gfa_file_name) as gfa_file:
        for gfa_line in gfa_file.readlines():
            # S       13      ctagtgggggacaaccaatcactaattgtgataAgatgtggtgtacacacactggtactggtcaggcaat  AS:37h; SP:0_70;        BP:0_61;        EV:1    FC:i:60 min:113 max:113 mean:113.0      AC:113;
            if gfa_line[0] != "S":
                continue
            print(f">{gfa_line.split()[1]}\n{gfa_line.split()[2]}")


if __name__ == '__main__':
    gfa_to_fa(sys.argv[1])
