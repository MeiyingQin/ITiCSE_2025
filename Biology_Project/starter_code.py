"""
DNA and Protein Mutation Analysis Tool
======================================

This Python module provides a functionalities for identifying, analyzing, and visualizing
mutations in DNA sequences. It covers mutations at both the DNA level (insertions, deletions,
and substitutions) and the amino acid level (missense, nonsense, no change). The tool can be
applied to compare wild-type and patient DNA sequences, with the assumption that the patient
sequence differs by a single base mutation.

Key Features:
-------------
1. **Sequence Validation (verify_dna_sequence):**
   - Validates whether a given sequence is a valid DNA (A, T, C, G) sequence.

2. **Base Distribution Calculation (get_dna_distribution):**
   - Computes the distribution of nucleotide bases (A, T, C, G) in a given DNA sequence.

3. **Open Reading Frame (ORF) Extraction (get_orf):**
   - Identifies all possible open reading frames (ORFs) from an RNA sequence.
   - Extracts and returns the longest ORF found in the sequence.

4. **DNA Transcription (transcribe_dna):**
   - Transcribes a DNA sequence (template sequence) into its corresponding RNA sequence.

5. **DNA Mutation Detection (identify_dna_mutation):**
   - Detects mutations at the DNA level by comparing the wild-type and patient sequences.
   - Classifies mutations into one of three categories: insertion, deletion, or substitution.
   - Provides details on the type of mutation, its position, and the base involved.

6. **Amino Acid Mutation Detection (identify_aa_mutation):**
   - Identifies mutations at the amino acid level by comparing the wild-type and patient RNA sequences.
   - Classifies mutations into missense (amino acid change) or nonsense (premature stop codon) mutations.
   - Provides information about the affected amino acid, its position, and the type of mutation.

7. **BRCA1 Mutation Check (is_brac_mutation):**
   - Confirms whether the patient's sequence differs from the wild type by exactly one base mutation.

8. **Mutation Visualization:**
   - **DNA-Level Visualization (visualize_dna_mutation):**
     - Visually represents DNA mutations (insertion, deletion, substitution) in a clear format, focusing
       on the mutation site and nearby bases.
   - **Amino Acid-Level Visualization (visualize_aa_mutation):**
     - Visually represents amino acid mutations, showing the wild-type and patient protein sequences and
       highlighting any changes.

9. **File Reading (read_sequence_from_file):**
    - Reads and retrieves DNA sequences from files, ensuring they are correctly formatted.

10. **Mutation Report (main):**
    - Performs the mutation analysis between the wild-type and patient DNA sequences.
    - Generates a detailed report, including:
      - The type of DNA mutation (if any) and its visualization.
      - The type of amino acid mutation (if any) and its visualization.
      - A comparison of the ORFs from the wild-type and patient sequences, highlighting any differences
        in start/stop codons and mutation locations within the ORF.

Workflow:
---------
1. The wild-type and patient DNA sequences are read from external files.
2. Sequences are validated to ensure they are valid DNA sequences.
3. A preliminary test was run to make sure the patient DNA is similar to the wild type
3. DNA sequences are transcribed into RNA.
4. The longest open reading frame (ORF) is extracted from the wild-type
5. The ORF on the patient RNA sequences is found based on the ORF of the wild-type.
5. DNA-level mutations (insertion, deletion, substitution) are identified and classified.
6. Amino acid-level mutations (missense, nonsense, etc.) are detected based on the DNA mutations.
7. Both DNA and amino acid mutations are visualized, focusing on relevant regions of the sequences.
8. A mutation report is generated, comparing the wild-type and patient sequences and
   explaining the potential effects of the mutation.

Assumptions:
------------
- Only a single base mutation is present between the wild-type and patient sequences.
- Mutations are classified as either insertion, deletion, or substitution at the DNA level.

Author: Amirhossein Nadiri
Date: 10/2/2024
Version: 1.0
"""
# pylint: disable=too-many-locals,too-many-arguments,too-many-function-args,too-many-branches
from enum import Enum
import doctest

# Enum for amino acid mutation types
class AAMutation(Enum):
    """
    AAMutation is an enumeration that represents the possible types of amino acid mutations
    that can occur due to a mutation in the DNA sequence.

    Types:
    - SILENT: No amino acid mutation.
    - MISSENSE: A mutation that changes only one of the amino acids to another amino acid.
    - NONSENSE: A mutation that only results in substitution of an amino acid to an early stop codon.
    - FRAMESHIFT: Multiple amino acid changes due to insertion/deletion, with possibility of an early stop.
    """
    SILENT = 0
    MISSENSE = 1
    NONSENSE = 2
    FRAMESHIFT = 3


class DNAMutation(Enum):
    """
    DNAMutation is an enumeration that represents different types of mutations
    that can occur in DNA sequences.

    Types:
    - NONE: No mutation is present, meaning the DNA sequences are identical.
    - INSERTION: A mutation where an extra base is added to the DNA sequence.
    - DELETION: A mutation where a base is removed from the DNA sequence.
    - SUBSTITUTION: A mutation where one base is replaced by another.
    - MULTI: Multiple mutations are present.
    """
    NONE = 0
    INSERTION = 1
    DELETION = 2
    SUBSTITUTION = 3
    MULTI = 4


codon_map = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '-', 'UAG': '-',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '-', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def verify_dna_sequence(sequence: str) -> bool:
    """
    Validate if a given string is a valid DNA sequence based on its type.

    Parameters:
    - sequence (str): The DNA sequence to validate.

    Returns:
    - bool: True if the sequence contains valid bases, False otherwise.

    Example:
    >>> verify_dna_sequence("ATCGTACG")
    True
    """
    # Please complete this function
    pass


def get_dna_distribution(sequence: str) -> dict:
    """
    Calculate the base distribution (A, G, C, T) in a given DNA sequence.

    Parameters:
    - sequence (str): The DNA sequence.

    Returns:
    - dict: A dictionary with base counts (A, G, C, T).

    Example:
    >>> get_dna_distribution("AGCTAGCT")
    {'A': 2, 'G': 2, 'C': 2, 'T': 2}
    """
    # Please complete this function
    pass


def get_orf(sequence: str) -> str:
    """
    Extract the longest open reading frame (ORF) from the RNA sequence.

    Parameters:
    - sequence (str): The RNA sequence.

    Returns:
    - str: The longest ORF found in the sequence.
    """
    # Please complete this function
    pass


def transcribe_dna(sequence: str) -> str:
    """
    Transcribe a given DNA sequence into RNA.

    Parameters:
    - sequence (str): The DNA sequence.

    Returns:
    - str: The transcribed RNA sequence.

    Example:
    >>> transcribe_dna("ATGC")
    'UACG'
    """
    # Please complete this function
    pass


def is_brac_mutation(wild_type: str, patient: str) -> bool:
    """
    Determine whether the patient sequence is a mutation of the wild-type DNA sequence,
    allowing for the possibility that the patient sequence may be incomplete (i.e.,
    some bases from its ending have been skipped).

    Parameters:
    - wild_type (str): The wild-type DNA sequence.
    - patient (str): The patient's DNA sequence.

    Returns:
    - bool: True if the ATCG distribution of the patient sequence is similar to the wild type. 
    Please note that the wild type sequence and the patient sequence may not have the same length. For the longer sequence, you take
    take the subsequence it (start from index 0) with the same length as the shorter sequence. Then compare the distribution of the two
    sequences. If the nucleotide count is different by less than 1 for all nucleotides. Then the two sequences are similar.
    
    Example:
    >>> is_brac_mutation("ATGCCGAC", "ATGGCGAC")
    True
    """
    # Please complete this function
    pass


def get_corresponding_orf(wild_type: str, patient: str) -> str:
    """
    Get the ORF in the patient sequence that corresponds to the ORF in the wild-type sequence.

    Assumption: There is only one mutation in the ORF, and the mutation is not at the start or
    stop codons of the wild-type ORF. The mutation may lead to an early stop.

    Parameters:
    - wild_type (str): The wild-type RNA sequence.
    - patient (str): The patient's RNA sequence.

    Returns:
    - str: The corresponding ORF in the patient sequence. If no valid ORF is found, return an empty string.
    """
    # Please complete this function
    pass


def identify_dna_mutation(wild_type: str, patient: str) -> (DNAMutation, int):
    """
    Identify the type and location of DNA mutation between a wild-type and patient sequence.

    Parameters:
    - wild_type (str): The wild-type DNA sequence.
    - patient (str): The patient's DNA sequence.

    Returns:
    - tuple: (DNAMutation type, mutation index).

    Example:
    >>> identify_dna_mutation("AACATACG", "AACCTACG")
    (<DNAMutation.SUBSTITUTION: 3>, 3)
    """
    # Please complete this function
    pass


def translate_rna(orf: str) -> str:
    """
    Translate a given RNA sequence (ORF) into its corresponding amino acid sequence.

    Parameters:
    - orf (str): An RNA sequence consisting of codons (e.g., 'AUGGCU...').

    Returns:
    - str: A string representing the amino acid sequence (e.g., 'MA...').
    """
    # Please complete this function
    pass


def identify_aa_mutation(
    wild_type: str,
    patient: str,
    dna_mutation: DNAMutation
) -> (list, AAMutation):
    """
    Identify all amino acid mutations and classify their types between the wild-type and patient protein sequences.

    Handles:
    - Silent: No amino acid mutation.
    - Missense: A mutation that changes only one of the amino acids to another amino acid.
    - Nonsense: A mutation that only results in substitution of an amino acid to an early stop codon.
    - FrameShift: Multiple amino acid changes due to insertion/deletion, with possibility of an early stop.

    Parameters:
    - wild_type (str): The wild-type protein sequence (without stop codons).
    - patient (str): The patient's protein sequence (without stop codons).
    - dna_mutation (DNAMutation): The type of DNA mutation (insertion, deletion, substitution).

    Returns:
    - tuple: (dictionary of mutation tuples (wild_aa, patient_aa), mutation classification as AAMutation).

    Example:
    >>> identify_aa_mutation("ML", "MI", DNAMutation.SUBSTITUTION)
    ({1: ('L', 'I')}, <AAMutation.MISSENSE: 1>)
    >>> identify_aa_mutation("MILIYILI", "MILI", DNAMutation.SUBSTITUTION)
    ({4: ('Y', '*'), 5: ('I', '*'), 6: ('L', '*'), 7: ('I', '*')}, <AAMutation.NONSENSE: 2>)
    """
    # Please complete this function
    pass    

def visualize_dna_mutation(
    wild_type: str,
    patient: str,
    mutation_type: DNAMutation,
    mutation_index: int,
    window: int = -1
) -> str:
    """
    Visualize the DNA mutation between the wild-type and patient sequences, adding ellipses (...)
    if parts of the sequence are not shown.

    Parameters:
    - wild_type (str): The wild-type DNA sequence.
    - patient (str): The patient's DNA sequence.
    - mutation_type (DNAMutation): The type of mutation.
    - mutation_index (int): The position of the mutation in the sequence.
    - window (int): The number of bases to show before and after the mutation (default: -1, no limit).

    Returns:
    - str: A formatted string visualizing the DNA mutation with a relevant portion of the sequence.
    
        If you want to generate something fancier than a string, you are welcome to!! If you want to import external modules for this,
    Please do not include it in the autograder as it may not be able to find the module and may raise errors. (You may remove this function in this case).
    """
    # Please complete this function
    pass


def visualize_aa_mutation(
    wild_type: str,
    patient: str,
    mutation_type: AAMutation,
    mutation_dict: dict,
    window: int = -1
) -> str:
    """
    Visualize the amino acid mutation between the wild-type and patient protein sequences.

    Parameters:
    - wild_type (str): Wild-type amino acid sequence (as a string of single-letter amino acids).
    - patient (str): Patient's amino acid sequence (as a string of single-letter amino acids).
    - mutation_type (AAMutation): The type of amino acid mutation (missense, nonsense, or frameshift).
    - mutation_dict (dict): A dictionary of mutation tuples containing (wild_aa, patient_aa).
    - window (int): The number of amino acids to show before and after the mutation (default: -1, no limit).

    Returns:
    - str: A formatted string representing the mutation visualization and its type.
    
    If you want to generate something fancier than a string, you are welcome to!! If you want to import external modules for this,
    Please do not include it in the autograder as it may not be able to find the module and may raise errors. (You may remove this function in this case).
    """
    # Please complete this function
    pass



def read_sequence_from_file(filename: str) -> str:
    """
    Read a DNA or RNA sequence from a file and return it.

    Parameters:
    - filename (str): The path to the file containing the sequence.

    Returns:
    - str: The sequence as a string if valid.
    """
    # Please complete this function
    pass


def generate_report(
    dna_mutation: DNAMutation,
    dna_visualization: str,
    aa_mutation: AAMutation,
    mutation_dict: dict,
    wild_type_protein: str,
    patient_protein: str,
    aa_visualization: str,
    wild_type_distribution: dict
) -> str:
    """
    Generate a detailed report of the mutation analysis at both the DNA and protein levels.

    Parameters:
    - dna_mutation: The type of mutation at the DNA level.
    - dna_visualization: Visualization of the DNA mutation.
    - aa_mutation: The type of mutation at the amino acid level.
    - wild_type_protein: The wild-type protein sequence.
    - patient_protein: The patient protein sequence.
    - aa_visualization: Visualization of the amino acid mutation.
    - wild_type_distribution: A dictionary representing the base counts in the wild-type DNA sequence.

    Returns:
    - str: The generated report as a string.
    """
    # Calculate the length difference and percentage change in the protein sequence
    protein_length_change = len(patient_protein) - len(wild_type_protein)
    protein_percentage_change = (
        (protein_length_change / len(wild_type_protein)) * 100 if len(wild_type_protein) > 0 else 0
    )

    # Initialize report as a list of strings to concatenate later.
    report = []

    # Include wild-type DNA nucleotide distribution in the report.
    report.append("Wild-Type DNA Nucleotide Distribution:")
    report.append(f"  Adenine (A): {wild_type_distribution['A']}")
    report.append(f"  Thymine (T): {wild_type_distribution['T']}")
    report.append(f"  Cytosine (C): {wild_type_distribution['C']}")
    report.append(f"  Guanine (G): {wild_type_distribution['G']}")

    # BRCA1 gene patient sequence section.
    report.append("\n-----------------------------------\n")
    report.append("BRCA1 Gene (Patient Sequence)")
    report.append("Changed from Wild-Type: " + ("Yes" if dna_mutation != DNAMutation.NONE else "No"))
    report.append(f"Type of DNA Mutation: {dna_mutation.name.title()}")

    # Print DNA alignment.
    report.append("DNA Sequence Alignment:")
    report.append(dna_visualization)

    # BRCA1 protein patient sequence section.
    report.append("-----------------------------------\n")
    report.append("BRCA1 Protein (Patient Sequence)")
    report.append("Change from Wild-Type: " + ("Yes" if aa_mutation != AAMutation.SILENT else "No"))
    report.append(f"Type of Amino Acid Mutation: {aa_mutation.name.title()}")

    # Report protein length absolute change.
    report.append(f"Change in Protein Length: {'Yes' if protein_length_change != 0 else 'No'} ({protein_length_change} amino acids)")

    # Determine if the ORF length has changed.
    if len(patient_protein) > len(wild_type_protein):
        report.append("The ORF in the patient sequence is longer.\n")
    elif len(patient_protein) < len(wild_type_protein):
        report.append("The ORF in the patient sequence is shorter.\n")

    # Report protein length percentage change.
    report.append(f"Percent Change in Protein Length: {abs(protein_percentage_change):.2f}%")

    # Length ratio.
    if len(wild_type_protein) > 0:
        length_ratio = len(patient_protein) / len(wild_type_protein)
        report.append(f"Protein Length Ratio (Patient/WT): {length_ratio:.2f}")

    # Protein consequence report.
    if aa_mutation == AAMutation.SILENT:
        report.append("Mutation:       No amino acid mutation (sequences are identical).")
    else:
        mutation_index, (wild_type_aa, patient_aa) = min(mutation_dict.items(), key=lambda x: x[0])
        if aa_mutation == AAMutation.MISSENSE:
            report.append(f"Mutation:       Missense mutation at position {mutation_index + 1}.")
            report.append( f"Amino acid '{wild_type_aa}' changed to '{patient_aa}'")
        elif aa_mutation == AAMutation.NONSENSE:
            report.append(f"Mutation:       Nonsense mutation at position {mutation_index + 1}.")
            report.append("A stop codon was introduced, leading to early termination of the protein.")
            removed_aas = [mutation[0] for mutation in mutation_dict.values()]
            removed_aas_string = ''.join(removed_aas)
            report.append(f"Amino acid(s) '{removed_aas_string}' were removed.")
        elif aa_mutation == AAMutation.FRAMESHIFT:
            report.append(f"Mutation:       Frameshift mutation starting at position {mutation_index + 1}.")
            report.append("Mismatched amino acids (index: WT -> Patient):")
            for mutation_index, (wt_aa, patient_aa) in mutation_dict.items():
                report.append(f"  - Position {mutation_index + 1}: {wt_aa} -> {patient_aa}")

    aa_difference_count = len(mutation_dict)
    aa_percentage_change = (aa_difference_count / len(wild_type_protein)) * 100
    report.append(f"Percent Different Amino Acids: {aa_percentage_change:.2f}%")

    # Print amino acid sequence alignment.
    report.append("Amino Acid Sequence Alignment:")
    report.append(aa_visualization)

    # Join the report list into a single string and return it.
    return "\n".join(report)


def main(wild_type_file: str, patient_file: str) -> str:
    """
    Perform the complete mutation analysis between wild-type and patient DNA sequences.

    Parameters:
    - wild_type_file (str): File path of the wild-type DNA sequence.
    - patient_file (str): File path of the patient's DNA sequence.

    Returns:
    - str: A report detailing the mutation analysis.
    """
    # Read the WT sequence from a text file.
    wild_type_dna = read_sequence_from_file(wild_type_file)

    # Verify that it is a DNA sequence as opposed to RNA or other
    # invalid content or content in inappropriate format etc.
    assert verify_dna_sequence(wild_type_dna), "Error: Invalid DNA sequences provided."

    # Find the distribution of the nucleotides.
    wild_type_distribution = get_dna_distribution(wild_type_dna)

    # Transcribe the DNA into an RNA sequence.
    wild_type_rna = transcribe_dna(wild_type_dna)

    # Find the open reading frame (ORF) of the sequence.
    wild_type_orf = get_orf(wild_type_rna)

    # Find the translated protein based on the ORF.
    wild_type_protein = translate_rna(wild_type_orf)

    # Read the patient sequence from a text file.
    patient_dna = read_sequence_from_file(patient_file)

    # Verify that it is a DNA sequence.
    assert verify_dna_sequence(patient_dna), "Error: Invalid patient DNA sequences provided."

    # Verify that that the patient sequence is a BRAC1 variant sequence.
    assert is_brac_mutation(wild_type_dna, patient_dna), "No single valid mutation detected (more than one base difference or no base difference)."

    # Transcribe the DNA into an RNA sequence.
    patient_rna = transcribe_dna(patient_dna)

    # On the mutation patient sequence, find the ORF subsequence that
    # corresponds to the ORF of the wild-type.
    patient_orf = get_corresponding_orf(wild_type_rna, patient_rna)

    # Find the type of DNA mutation.
    dna_mutation, mutation_index = identify_dna_mutation(wild_type_dna, patient_dna)

    # Translate the patient ORF into protein.
    patient_protein = translate_rna(patient_orf)

    # Identify the consequences on the patient protein.
    mutation_dict, mutation_type = identify_aa_mutation(wild_type_protein, patient_protein, dna_mutation)

    # Generate a visual aligning the WT DNA with that of the patient.
    dna_visualization = visualize_dna_mutation(wild_type_dna, patient_dna, dna_mutation, mutation_index, 10)

    # Generate a visual aligning the WT protein sequence with that of the patient.
    aa_visualization = visualize_aa_mutation(wild_type_protein, patient_protein, mutation_type, mutation_dict, 10)

    # Report Generation
    return generate_report(
        dna_mutation,
        dna_visualization,
        mutation_type,
        mutation_dict,
        wild_type_protein,
        patient_protein,
        aa_visualization,
        wild_type_distribution
    )


# Example usage:
if __name__ == "__main__":
    doctest.testmod()
    # Please put the sequences (e.g., wild_type.txt, patient_1.txt) in the same directory/folder as this Python script
    ANALYSIS_REPORT = main('wild_type.txt', 'patient_1.txt') # You may change the file name to see the results of other patients
    print(ANALYSIS_REPORT)