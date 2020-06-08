import sys
"""
Given a fasta file containing n sequences having all the same length
Returns the pi value https://en.wikipedia.org/wiki/Nucleotide_diversity


1/ stores the n sequences
2/ computes for all i<j in [1,n] the pi_i,j value
3/ computes and prints the pi value
"""


def stores_fasta_sequences(fasta_file_name):
    sequences=[]
    with  open(fasta_file_name) as fasta_file:
        while True: 
            line = fasta_file.readline()
            if not line: break                  # end of the file
            line = line.strip()
            if line[0] == ">": continue         # comment
            sequences.append(line)
    if len(sequences) == 0: return None
    l = len(sequences[0])
    for sequence in sequences:
        if len(sequence) != l: 
            sys.stderr.write(f"File {fasta_file_name} contains sequences of distinct sizes, unable to compute pi\n")
            sys.exit(1)
    return sequences
    

def main():
    if len(sys.argv) !=2:
        sys.stderr.write(f"Usage: {sys.argv[0]} fasta_file_name\n")
        sys.stderr.write("Given a fasta file containing n sequences having all the same length\n")
        sys.stderr.write("Returns the pi value https://en.wikipedia.org/wiki/Nucleotide_diversity\n")
        sys.stderr.write(" 1/ stores the n sequences\n")
        sys.stderr.write(" 2/ computes for all i<j in [1,n] the pi_i,j value\n")
        sys.stderr.write(" 3/ computes and prints the pi value\n")
    
    sequences = stores_fasta_sequences(sys.argv[1])
    if not sequences :
        sys.stderr.write(f"File {sys.argv[1]} contains no sequence\n")
    print (sequences)

if __name__ == '__main__':
    main()

