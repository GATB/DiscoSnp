import sys
# import pickle
import getopt

    
comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}


def rev_comp(seq):
    return (''.join(comp[c] for c in seq))[::-1]
    


def get_huigher_only(seq):
    res=""
    for l in seq:
        if l>='A' and l<='Z': 
            res+=l
    return res
    
def lower_central(seq):
    m=len(seq)//2
    return seq[:m]+seq[m].lower()+seq[m+1:]

def index_fasta(fasta_file_name,id2seq = {}):
    fasta_file = open(fasta_file_name)
    while True:
        com= fasta_file.readline() #>SNP_higher_path_911068|P_1:30_C/T|high|nb_pol_1|left_unitig_length_9|right_unitig_length_25|C1_115|C2_0|Q1_69|Q2_0|G1_0/0:9,350,2304|G2_1/1:84,16,4|rank_1
        if not com: break
        seq= get_huigher_only(fasta_file.readline())
        idsnp=com.split('|')[0].split('_')[3]+com.split('|')[0].split('_')[1][0] # 911068h
        id2seq[idsnp]=lower_central(seq)
        
        
    fasta_file.close()
    return id2seq
    

def print_consensus_fact(fact_file_name,id2seq):
    fact_file = open(fact_file_name)
    fact_id=-1
    previous_sequence_length=0
    while True:
        line = fact_file.readline()
        if not line: break
        if not line.startswith("fact"): continue
        this_fact_id = line.split(',')[1] #fact(cc4,0,1,1000020,m,h,0).
        if this_fact_id != fact_id:
            print("\n>FACT_ID_"+str(this_fact_id))
            fact_id=this_fact_id
            previous_sequence_length=0
            
        snp_id=line.split(',')[3]+line.split(',')[5]
        rc = False
        sign=""
        if line.split(',')[4] == 'n':
        # if snp_id[0]=='-':
            rc=True
            sign='-'
        if snp_id not in id2seq: # IN UNCOHERENT
            continue
        else: seq = id2seq[snp_id]
        
        
        if rc: 
            seq = rev_comp(seq)
        shift=previous_sequence_length-int(line.split(',')[-1].split(')')[0])
        previous_sequence_length=len(seq)
        if shift<0:
            for i in range(-shift): print('N',end='')
            print (seq,end='')
        else:
            print (seq[shift:],end='')
            
        
        
        
        
def print_fact(fact_file_name,id2seq):
    fact_file = open(fact_file_name)
    fact_id=-1
    cumulated_shift=0
    while True:
        line = fact_file.readline()
        if not line: break
        if not line.startswith("fact"): continue
        this_fact_id = line.split(',')[1] #fact(cc4,0,1,1000020,m,h,0).
        if this_fact_id != fact_id:
            print("\n FACT ID",this_fact_id)
            fact_id=this_fact_id
            cumulated_shift=0
        snp_id=line.split(',')[3]+line.split(',')[5]
        rc = False
        sign=""
        if line.split(',')[4] == 'n':
        # if snp_id[0]=='-':
            rc=True
            sign='-'
        if snp_id not in id2seq: # IN UNCOHERENT
            seq=None
            cumulated_shift+=int(line.split(',')[-1].split(')')[0])
            continue
        else: seq = id2seq[snp_id]
        
        
        if rc: 
            seq = rev_comp(seq)
        # alt seq: 
        if snp_id[-1]=='h':
            alt_seq = id2seq[snp_id[:-1]+'l']
        else:
            alt_seq = id2seq[snp_id[:-1]+'h']
        if rc: alt_seq = rev_comp(alt_seq)

            
        seq+=" ("+alt_seq[len(alt_seq)//2]+")"+" "+sign+snp_id
        cumulated_shift+=int(line.split(',')[-1].split(')')[0])
        shift = cumulated_shift         
        for i in range(shift):
            print(" ",end='')
        print (seq)
    
def usage():
    usage= """
    #########################################
    facts2aligned_sequences.py
    #########################################
    
    -h --help : print this message
    --coherent_file:            <file>.fa:  coherent fa file from discoSnp    [Mandatory]
    --uncoherent_file:          <file>.fa:  uncoherent fa file from discoSnp  [Optional]
    --fact_file:                <file>.txt: file generated from format_phased_variants_for_haplotyping.py
    --consensus:                Boolean:    It set: print a consensus per fact instead of stacked sequences
    """
    print(usage,file=sys.stderr)
def main():
    uncoherent_fa_file=None
    print_consensus=False
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:u:f:C",["help","coherent_file=","uncoherent_file=","fact_file=","consensus"])
        if not opts:
            usage()
            sys.exit(2)
    except getopt.GetoptError as e:
        print(e)
        usage()
        sys.exit(2)
    for opt, arg in opts : 
        if opt in ("-h", "--help"):
            usage()
            sys.exit(2)
        elif opt in ("-c","--coherent_file"):
            coherent_fa_file_name = arg
        elif opt in ("-u","--uncoherent_file"):
            uncoherent_fa_file_name = arg
        elif opt in ("-f","--fact_file"):
            fact_file_name = arg
        elif opt in ("-C","--consensus"):
            print_consensus=True
        else:
            print("Unkwnown option {} ".format(opt),file=sys.stderr)
            usage()
            sys.exit(2)
            
            
    id2seq=index_fasta(coherent_fa_file_name)
    if uncoherent_fa_file: id2seq=index_fasta(uncoherent_fa_file_name,id2seq)
# pickle.dump(id2seq, open("id2seq","wb"))
    
    # id2seq=pickle.load(open("id2seq","rb"))
    
    
    
        #
    # if sys.argv[1][-1]=='a':
    #     id2seq=index_fasta(sys.argv[1])
    #     pickle.dump(id2seq, open("id2seq","wb"))
    # else:
    #     id2seq=pickle.load(open("id2seq","rb"))
    if print_consensus:
        print_consensus_fact(fact_file_name,id2seq)
    else:
        print_fact(fact_file_name,id2seq)
if __name__ == '__main__':
    main()