# =================
# author: Sara A. Carioscia, Biology Dept., Johns Hopkins University
# email: scarios1@jhu.edu
# last update: June 14, 2024
# aim: This script reads a map file containing genomic position and chromosome (two columns) for each site in a vcf,
#       divides the data into a specified number of chunks, and writes the chunk information to an output file.
#       Each chunk includes the chromosome and position range for the lines in that chunk.
# =================

# Usage:
#     python3 create_chunks.py <num_chunks> <input_mapfile> <outfilename>

#     num_chunks       : Number of chunks to divide the input file into.
#     input_mapfile    : Path to the input file containing genomic data.
#     outfilename      : Path to the output file where chunk information will be written.

import sys
import math

def create_chunks(num_chunks, input_file, output_file):
	# Read the lines from the input file
	with open(input_file, 'r') as f:
		lines = f.readlines()

	# Calculate the number of lines per chunk
	lines_per_chunk = math.ceil(len(lines) / num_chunks)

	# Create chunks and write to the output file
	with open(output_file, 'w') as out_file:
		# Iterate over each chunk
		for i in range(num_chunks):
			# Calculate the start and end positions for the current chunk
			start_idx = i * lines_per_chunk
			end_idx = min(start_idx + lines_per_chunk, len(lines))

			# Extract the lines for the current chunk
			chunk_lines = lines[start_idx:end_idx]

			# Extract the chromosome and positions from the chunk lines
			chromosome = chunk_lines[0].split()[0]
			start_pos = chunk_lines[0].split()[1]
			end_pos = chunk_lines[-1].split()[1]

			# Format the chunk information
			chunk_info = f"{chromosome}:{start_pos}-{end_pos}\n"

			# Write the chunk information to the output file
			out_file.write(chunk_info)


if __name__ == "__main__":
	# Get the number of chunks and input file from command line arguments
	if len(sys.argv) != 4:
		print("Usage: python script.py <num_chunks> <input_mapfile> <outfilename>")
		sys.exit(1)

	num_chunks = int(sys.argv[1])
	input_file = sys.argv[2]
	output_file = sys.argv[3]

	# Create the chunks
	create_chunks(num_chunks, input_file, output_file)
