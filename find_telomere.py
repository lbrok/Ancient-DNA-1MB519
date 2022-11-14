def find_telomere():
    # Reading fasta file to string
    fasta_file = open(r"C:/Users/46722/Documents/Tillämpad_bioinformatik/Ancient-DNA-1MB519/chr13.fa", "r")
    seq = fasta_file.read()
    fasta_file.close()

    # Removing first line of fasta sequence (identifier)
    seq = seq.split('\n',2)[-1]

    # Finding first position that is not N
    A_pos = seq.find('A')
    G_pos = seq.find('G')
    C_pos = seq.find('C')
    T_pos = seq.find('T')
    first_pos = min(A_pos, G_pos, C_pos, T_pos)

    return first_pos

if __name__ == '__main__':
    print(f"The last chromosome position of the telomere is: {find_telomere()-1} bp")

