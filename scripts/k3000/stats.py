import matplotlib.pyplot as plt
import sys


def get_nb_allele_distribution(gfa_file_name):
    gfa_file=open(gfa_file_name)
    sizes = []
    for line in gfa_file:
        if line[0]=='S':#S  2       34156l;-11363l;13698l;-26143h;10014l;   RC:i:144
            l=len(line.split()[2].split(';')[:-1])
            sizes.append(l)
    gfa_file.close()
    return sizes

def get_sequence_size_distribution(gfa_file_name):
    gfa_file=open(gfa_file_name)
    sizes = []
    for line in gfa_file:
        if line[0]=='S':#S       0       actgcaACAGCTGTTGAAAAGCCGGAATGTACTCTTCATTGCAAACATTTCAGGGATGAAGTGAAGAtgaattgCGACGTAGTATCCACACCAAGCCGGCGTTATCCGGTGAGGCGCAATGTTGCGGGGGCttt  RC:i:11
            sizes.append(len(line.split()[2]))
    gfa_file.close()
    return sizes
    
def plot_violin(sequence_sizes, nb_alleles):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    axes[0].violinplot(sequence_sizes)
    axes[0].set_title('Sequence size distribution')
    axes[1].violinplot(nb_alleles)
    axes[1].set_title('Nb alleles per sequence')
    plt.savefig('distributions.png')
    
    

def main():
    '''
    Stats from a gfa file
    '''
    sequence_sizes=get_sequence_size_distribution(sys.argv[2]) #for each snp id: sequences[snp_id]=[left_unitig_len, right_unitig_len, upperseq, lowerseq] 
    nb_alleles=get_nb_allele_distribution(sys.argv[1])
    plot_violin(sequence_sizes, nb_alleles)
    
    



if __name__ == "__main__":
     main()
