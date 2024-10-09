import pysam
import csv
import sys

bam = sys.argv[1]
id = sys.argv[2]

output_file = f"{id}_depth_analysis.csv"

def analyze_bam(bam, output_file):
    depth_data = []
    print("Analyzing BAM file for depth information...")

    # Open BAM file and iterate through positions
    with pysam.AlignmentFile(bam, "rb") as bam_file:
        for pileupcolumn in bam_file.pileup():
            position = pileupcolumn.pos + 1  # 1-based position
            base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0}

            # Count base occurrences at each position
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.query_position is not None:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] = base_counts.get(base, 0) + 1
                else:
                    # Count reads that don't align to a specific base as 'N'
                    base_counts['N'] += 1

            # Use the sum of base-specific depths as the total depth
            total_depth = sum(base_counts.values())

            # Append data for CSV
            depth_data.append([position, total_depth] + list(base_counts.values()))

    # Write results to CSV
    with open(output_file, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Position', 'Total Depth', 'A Depth', 'T Depth', 'C Depth', 'G Depth', 'N Depth'])
        writer.writerows(depth_data)

    print(f"BAM depth analysis completed and saved to {output_file}")

# Run the function
analyze_bam(bam, output_file)
