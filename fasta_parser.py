def parse_fasta(file_path):
    """Parse a FASTA file and return a dictionary of sequence IDs and sequences."""
    sequences = {}
    with open(file_path, 'r') as file:
        sequence_id = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = "".join(sequence)
                sequence_id = line[1:]  # Remove '>' to get the sequence ID
                sequence = []
            else:
                sequence.append(line)
        if sequence_id:
            sequences[sequence_id] = "".join(sequence)  # Add the last sequence
    return sequences


def calculate_gc_content(sequence):
    """Calculate GC content as a percentage."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0


def main():
    file_path = input("Enter the path to the FASTA file: ")
    sequences = parse_fasta(file_path)
    for seq_id, sequence in sequences.items():
        gc_content = calculate_gc_content(sequence)
        print(f"{seq_id}: Length={len(sequence)}, GC%={gc_content:.2f}")


if __name__ == "__main__":
    main()
