# DNA-Restriction-Enzyme-Tool
A Python tool that analyzes nucleotide sequences from a FASTA file and identifies restriction enzyme cutting sites, generating a detailed report of resulting fragments.

## Features
- Reads a FASTA file and a list of restriction enzymes from command-line arguments.
- Identifies cutting sites for each enzyme in the sequence.
- Generates a detailed report including fragment lengths and formatted sequences.

## Usage
Run the script with the following command-line arguments:
```bash
python3 restriction_enzyme_analyzer.py <nucleotide_sequence.fasta> <restriction_enzymes.txt>
```
## Input Files

*   **FASTA File**: Contains a single nucleotide sequence (e.g., `input_sequence.fasta`).
*   **Enzyme File**: Lists restriction enzymes and their recognition sequences with a `^` (or `%`) indicating the cleavage site (e.g., `enzyme_list.txt`). The format should be one enzyme per line: `EnzymeName;Recognition^Sequence`.

## Output

The script prints a report to the terminal, including:

*   The names of the input files used.
*   The sequence header from the FASTA file.
*   The length of the nucleotide sequence.
*   For each enzyme:
    *   The name and recognition sequence.
    *   The number of cutting sites found.
    *   If sites are found: the number of fragments and their sequences, formatted in lines of 60 bases grouped into 10-base blocks with positional numbering.
    *   If no sites are found: a message stating this.

## Example

Example command:
```bash
python3 restriction_enzyme_analyzer.py input_sequence.fasta enzyme_list.txt
```
## Requirements

Python 3.x
No external dependencies beyond the standard library.
