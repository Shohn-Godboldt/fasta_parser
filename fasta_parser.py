import os
import matplotlib.pyplot as plt


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


def count_ambiguous_bases(sequence):
    """Count ambiguous bases (e.g., 'N') in the sequence."""
    return sequence.count('N')


def plot_gc_content(gc_values):
    """Plot the distribution of GC content."""
    plt.hist(gc_values, bins=10, color='blue', edgecolor='black')
    plt.title("GC Content Distribution")
    plt.xlabel("GC%")
    plt.ylabel("Frequency")
    plt.show()


def main():
    file_path = input("Enter the path to the FASTA file: ")
    
    # Error handling for file
    if not os.path.exists(file_path):
        print("Error: File not found.")
        return
    if not file_path.endswith(".fasta"):
        print("Error: File is not in FASTA format.")
        return
    
    sequences = parse_fasta(file_path)
    summary = []
    gc_values = []

    # Analyze sequences
    for seq_id, sequence in sequences.items():
        gc_content = calculate_gc_content(sequence)
        ambiguous_count = count_ambiguous_bases(sequence)
        gc_values.append(gc_content)
        summary.append(f"{seq_id}: Length={len(sequence)}, GC%={gc_content:.2f}, Ambiguous={ambiguous_count}")
    
    # Print results to terminal
    for line in summary:
        print(line)
    
    # Write results to a file
    with open("output.txt", "w") as output_file:
        output_file.write("\n".join(summary))
    
    print("\nAnalysis complete. Results saved to 'output.txt'.")
    
    # Plot GC content
    plot_gc_content(gc_values)


if __name__ == "__main__":
    main()

