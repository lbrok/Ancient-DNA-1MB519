from Bio.SeqIO.FastaIO import SimpleFastaParser

def find_telomere():
    # Reading fasta file to string
    with open("D:/Data/LoxAfr4_DQ188829.fa", encoding="ascii") as handle:
        for seq in SimpleFastaParser(handle):

            # Finding first position that is not N
            A_pos = seq[1].find('A')
            G_pos = seq[1].find('G')
            C_pos = seq[1].find('C')
            T_pos = seq[1].find('T')
            first_pos = min(A_pos, G_pos, C_pos, T_pos)
            print(f"The last chromosome position of the telomere in {seq[0]} is: {first_pos} bp")
            #return [seq[0],first_pos]

if __name__ == '__main__':
    result = find_telomere()
    #print(f"The last chromosome position of the telomere in {result[0]} is: {result[1]-1} bp")

