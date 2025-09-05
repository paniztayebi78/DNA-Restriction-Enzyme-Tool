import sys  # Imports sys module to handle command-line arguments and exit the program if arguments are missing
import re   # Imports the re module for regular expressions, which helps locate cutting sites within the sequence

def read_fasta(nucleotide_seq):
    """Reads the FASTA file to extract the header and nucleotide sequence.

    Args:
        nucleotide_seq (str): The file path of the FASTA file containing a single nucleotide sequence.

    Returns:
        tuple: The header of the sequence (str) and the nucleotide sequence (str) with all lines concatenated.
        
    This function opens the FASTA file, reads the header (first line), and concatenates the nucleotide sequence,
    removing newline characters for ease of further processing.
    """
    with open(nucleotide_seq, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip()[1:]  # Skip the initial '>' character in the FASTA header to get the sequence name
        # Concatenate remaining lines as a single string without newline characters
        sequence = ''.join(line.strip() for line in lines[1:])
    return header, sequence


def read_enzymes(restriction_enzymes):
    """Reads the enzyme file and extracts enzyme names and their recognition sequences with cutting sites.

    Args:
        restriction_enzymes (str): The file path of the enzyme file containing enzyme names and sequences.

    Returns:
        list: A list of tuples where each tuple contains an enzyme name and its recognition sequence (str).
        
    The enzyme file is expected to have one enzyme per line in the format "EnzymeName;RecognitionSequence^".
    This function maps each enzyme to its cutting pattern.
    """
    enzymes = []
    with open(restriction_enzymes, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines to avoid errors
                # Each line is split by the ';' character, separating enzyme name and sequence
                name, sequence = line.strip().split(';')
                enzymes.append((name, sequence))  # Append as a tuple (enzyme name, recognition sequence)
    return enzymes

def find_cuts(sequence, enzyme_sequence):
    """Finds the cutting positions within the nucleotide sequence for a specific enzyme.

    Args:
        sequence (str): The full nucleotide sequence where cuts will be found.
        enzyme_sequence (str): The recognition sequence of the enzyme, with '%' indicating the cutting position.

    Returns:
        list: A list of integer positions where the enzyme will cut the sequence.
        
    This function locates all positions in the sequence where the enzyme cuts, replacing '%' with '' in the pattern.
    It uses regular expressions to locate all matches of the recognition sequence within the nucleotide sequence.
    """
    cut_positions = []
    # Remove '%' to create a search pattern that matches the recognition sequence
    cut_pattern = enzyme_sequence.replace('%', '')
    cut_offset = len(cut_pattern)  # Offset to place cut right after the recognition sequence

    # Use regex to find matches of the cut pattern in the full nucleotide sequence
    for match in re.finditer(cut_pattern, sequence):
        cut_positions.append(match.start() + cut_offset)  # Append exact cut position
    return cut_positions

def generate_report(nucleotide_seq, restriction_enzymes, header, sequence, enzymes):
    """Generates and formats a report of enzyme cutting analysis.

    Args:
        nucleotide_seq (str): Name of the input FASTA file.
        restriction_enzymes (str): Name of the input enzyme file.
        header (str): Header of the nucleotide sequence from the FASTA file.
        sequence (str): The nucleotide sequence.
        enzymes (list): List of enzyme names and recognition sequences.

    Returns:
        str: The full formatted report as a single string.
        
    This function organizes the report to match the specified format, including details of each enzyme,
    where it cuts, fragment lengths, and formatted nucleotide fragments.
    """
    sequence_length = len(sequence)
    report = []

    # Initial report details with sequence file name, enzyme file name, and basic sequence info
    report.append(f"Restriction enzyme analysis of sequence from file {nucleotide_seq}.")
    report.append(f"Cutting with enzymes found in file {restriction_enzymes}.")
    report.append("-" * 63)
    report.append(f"Sequence name:  {header}")
    report.append(f"Sequence is {sequence_length} bases long.")
    report.append("-" * 63)
    
    # For each enzyme, find the cutting positions and generate fragment information
    for enzyme_name, enzyme_sequence in enzymes:
        cut_positions = find_cuts(sequence, enzyme_sequence)
        cut_count = len(cut_positions)
        
        report.append(f"There are {cut_count} cutting sites for {enzyme_name}, cutting at {enzyme_sequence}")
        
        if cut_count > 0:
            fragments = []
            start = 0
            # Create sequence fragments by slicing at each cutting position
            for cut in cut_positions:
                fragments.append(sequence[start:cut])  # Append each fragment from start to the cut position
                start = cut  # Update start to the current cut position for the next fragment
            fragments.append(sequence[start:])  # Add the final fragment after the last cut

            report.append(f"There are {len(fragments)} fragments:")

            # Format each fragment for output with positions and lengths
            position = 1  # Starting base position for display
            for fragment in fragments:
                length = len(fragment)
                # Create formatted lines for each fragment, splitting into lines of 60 bases grouped by 10
                fragment_lines = [f"{position:3} " + ' '.join([fragment[i:i+10] for i in range(0, len(fragment), 10)])
                                  for position in range(0, len(fragment), 10)]
                report.append(f"Length- {length}")
                report.extend(fragment_lines)
                position += length  # Update the position for the next fragment
        else:
            report.append(f"There are no sites for {enzyme_name}.")
        
        report.append("-" * 63)
    
    return "\n".join(report)  # Join all report lines into a single formatted string

def main():
    """Main function to process input files and display the restriction enzyme analysis report."""
    # Check for correct number of command-line arguments (expecting 2)
    if len(sys.argv) != 3:
        print("Usage: python3 Assignment2.py nucleotide_seq.fasta restriction_enzymes.txt")
        sys.exit(1)  # Exit if the required arguments are missing

    # Read filenames for FASTA and enzyme files from the command-line arguments
    nucleotide_seq = sys.argv[1]
    restriction_enzymes = sys.argv[2]

    # Step 1: Read the FASTA file to retrieve the header and sequence
    header, sequence = read_fasta(nucleotide_seq)

    # Step 2: Read the enzymes file to retrieve enzyme names and cutting patterns
    enzymes = read_enzymes(restriction_enzymes)

    # Step 3: Generate the report based on cutting sites and fragment information
    report = generate_report(nucleotide_seq, restriction_enzymes, header, sequence, enzymes)

    # Step 4: Print the final report to the console
    print(report)

if __name__ == "__main__":
    main()  # Execute the main function when the script is run directly
