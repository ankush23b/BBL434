import sys
from collections import Counter
import matplotlib.pyplot as plt

def read_fasta(filename):
    with open(filename, 'r') as f:
        # Skip header and join sequence lines
        return "".join(line.strip() for line in f if not line.startswith(">"))

def get_kmer_enrichment(sequence, k, window_size, step):
    enrichment_data = []
    max_count = 0
    best_loc = 0
    best_kmer = ""

    # Sliding window logic
    for i in range(0, len(sequence) - window_size + 1, step):
        window = sequence[i : i + window_size]
        
        # Count all kmers in the current window
        kmers = [window[j : j + k] for j in range(len(window) - k + 1)]
        counts = Counter(kmers)
        
        if not counts:
            continue
            
        # Find the most frequent k-mer in this specific window
        top_kmer, count = counts.most_common(1)[0]
        enrichment_data.append((i, count))
        
        # Track the global maximum for the final report
        if count > max_count:
            max_count = count
            best_loc = i
            best_kmer = top_kmer
            
    return enrichment_data, best_kmer, max_count, best_loc

def plot_enrichment(data, k, output_file):
    positions, counts = zip(*data)
    plt.figure(figsize=(12, 6))
    plt.plot(positions, counts, color='blue', linewidth=1)
    plt.title(f'K-mer Enrichment (k={k}) along Genome')
    plt.xlabel('Genomic Position (bp)')
    plt.ylabel('Highest K-mer Count in Window')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py genomic.fa")
        sys.exit(1)

    # Parameters
    K = 8
    WINDOW = 5000
    STEP = 500
    
    fasta_file = sys.argv[1]
    sequence = read_fasta(fasta_file)
    
    print(f"Processing sequence of length: {len(sequence)} bp")
    
    data, kmer, count, loc = get_kmer_enrichment(sequence, K, WINDOW, STEP)
    
    print("-" * 30)
    print(f"MOST ENRICHED LOCATION FOUND")
    print(f"Location (start): {loc} bp")
    print(f"K-mer sequence:   {kmer}")
    print(f"Occurrence count: {count}")
    print("-" * 30)
    
    plot_enrichment(data, K, "script1_plot1.png")

if __name__ == "__main__":
    main()
