import os
import matplotlib.pyplot as plt
import plotly.express as px
import pandas as pd


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


def parse_fastq(file_path):
    """Parse a FASTQ file and return sequences with their quality scores."""
    sequences = {}
    qualities = {}
    with open(file_path, 'r') as file:
        while True:
            seq_id = file.readline().strip()  # Sequence ID
            if not seq_id:
                break  # End of file
            sequence = file.readline().strip()  # Sequence
            file.readline()  # Plus line (skipped)
            quality = file.readline().strip()  # Quality scores
            sequences[seq_id[1:]] = sequence
            qualities[seq_id[1:]] = quality
    return sequences, qualities


def calculate_gc_content(sequence):
    """Calculate GC content as a percentage."""
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100 if sequence else 0


def calculate_avg_quality(quality):
    """Calculate average quality score from the quality string."""
    return sum(ord(char) - 33 for char in quality) / len(quality) if quality else 0


def plot_gc_content(gc_values):
    """Create a histogram of GC content."""
    plt.hist(gc_values, bins=10, color='blue', edgecolor='black')
    plt.title("GC Content Distribution")
    plt.xlabel("GC%")
    plt.ylabel("Frequency")
    plt.show()


def interactive_gc_plot(gc_values):
    """Create an interactive GC content plot."""
    df = pd.DataFrame({"GC Content": gc_values})
    fig = px.histogram(df, x="GC Content", nbins=10, title="GC Content Distribution",
                       labels={"GC Content": "GC%"},
                       template="plotly_dark")
    fig.show()


def main():
    file_path = input("Enter the path to the FASTA/FASTQ file: ")
    
    # Error handling for file existence and format
    if not os.path.exists(file_path):
        print("Error: File not found.")
        return
    if not file_path.endswith((".fasta", ".fastq")):
        print("Error: Unsupported file format. Only FASTA and FASTQ are supported.")
        return

    # Handle FASTA files
    if file_path.endswith(".fasta"):
        sequences = parse_fasta(file_path)
        gc_values = []
        summary = []

        for seq_id, sequence in sequences.items():
            gc_content = calculate_gc_content(sequence)
            gc_values.append(gc_content)
            summary.append(f"{seq_id}: Length={len(sequence)}, GC%={gc_content:.2f}")
        
        # Print and save results
        for line in summary:
            print(line)
        with open("output.txt", "w") as output_file:
            output_file.write("\n".join(summary))
        print("\nAnalysis complete. Results saved to 'output.txt'.")

        # Visualize GC content
        interactive_gc_plot(gc_values)

    # Handle FASTQ files
    elif file_path.endswith(".fastq"):
        sequences, qualities = parse_fastq(file_path)
        gc_values = []
        summary = []

        for seq_id, sequence in sequences.items():
            gc_content = calculate_gc_content(sequence)
            avg_quality = calculate_avg_quality(qualities[seq_id])
            gc_values.append(gc_content)
            summary.append(f"{seq_id}: Length={len(sequence)}, GC%={gc_content:.2f}, Avg Quality={avg_quality:.2f}")
        
        # Print and save results
        for line in summary:
            print(line)
        with open("output.txt", "w") as output_file:
            output_file.write("\n".join(summary))
        print("\nAnalysis complete. Results saved to 'output.txt'.")

        # Visualize GC content
        interactive_gc_plot(gc_values)


if __name__ == "__main__":
    main()
