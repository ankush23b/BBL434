import sys
import matplotlib.pyplot as plt

def read_fasta(filename):
    with open(filename, 'r') as f:
        return "".join(line.strip() for line in f if not line.startswith(">")).upper()

def calculate_gc_skew(sequence, window_size, step):
    positions = []
    cumulative_skews = []
    current_cumulative = 0
    
    # Track min/max for enrichment analysis
    min_skew_val = 0
    ori_index = 0

    for i in range(0, len(sequence) - window_size + 1, step):
        window = sequence[i : i + window_size]
        g = window.count('G')
        c = window.count('C')
        
        # Calculate skew for the window
        if (g + c) == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)
        
        current_cumulative += skew
        positions.append(i)
        cumulative_skews.append(current_cumulative)
        
        # The minimum of the cumulative plot is the likely OriC
        if current_cumulative < min_skew_val:
            min_skew_val = current_cumulative
            ori_index = i

    return positions, cumulative_skews, ori_index

def main():
    if len(sys.argv) < 2:
        print("Usage: python script3.py genomic.fa")
        sys.exit(1)

    # Parameters
    WINDOW = 5000
    STEP = 500
    
    sequence = read_fasta(sys.argv[1])
    positions, skews, ori_pos = calculate_gc_skew(sequence, WINDOW, STEP)
    
    print("-" * 35)
    print(f"GC SKEW ANALYSIS COMPLETE")
    print(f"Predicted OriC (Cumulative Min): {ori_pos} bp")
    print("-" * 35)

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.plot(positions, skews, color='green', linewidth=1.5)
    plt.axvline(x=ori_pos, color='red', linestyle='--', label=f'Predicted OriC: {ori_pos}')
    
    plt.title('Cumulative GC Skew Plot')
    plt.xlabel('Genome Position (bp)')
    plt.ylabel('Cumulative Skew')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig("script3_plot1.png")
    print("Plot saved as script3_plot1.png")

if __name__ == "__main__":
    main()
