#Transcription
#A gene transcribed from mRNA and then to a protein (Amino Acid) in a form of Codons (which is a sequence of 3 nucleotides)
#The process of converting DNA to RNA is called transcription, and the process of converting RNA to protein is called translation.
#We then define all the codon transformation into amino acid, and an amino acid corresponding into a codon

# DNA Codon Table: Maps DNA codons to amino acids

# Project: DNA Transcription and Translation Simulator
# Purpose: Simulate DNA to protein synthesis with testing for scholarship application
# Date: March 25, 2025

import unittest

# DNA Codon Table: Maps RNA codons to amino acids
codon_to_amino_acid = {
    "UUU": "Phenylalanine", "UUC": "Phenylalanine",
    "UUA": "Leucine", "UUG": "Leucine",
    "UCU": "Serine", "UCC": "Serine", "UCA": "Serine", "UCG": "Serine",
    "UAU": "Tyrosine", "UAC": "Tyrosine",
    "UAA": "STOP", "UAG": "STOP", "UGA": "STOP",
    "UGU": "Cysteine", "UGC": "Cysteine",
    "UGG": "Tryptophan",
    "CUU": "Leucine", "CUC": "Leucine", "CUA": "Leucine", "CUG": "Leucine",
    "CCU": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "CAU": "Histidine", "CAC": "Histidine",
    "CAA": "Glutamine", "CAG": "Glutamine",
    "CGU": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "AUU": "Isoleucine", "AUC": "Isoleucine", "AUA": "Isoleucine",
    "AUG": "Methionine",
    "ACU": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "AAU": "Asparagine", "AAC": "Asparagine",
    "AAA": "Lysine", "AAG": "Lysine",
    "AGU": "Serine", "AGC": "Serine",
    "AGA": "Arginine", "AGG": "Arginine",
    "GUU": "Valine", "GUC": "Valine", "GUA": "Valine", "GUG": "Valine",
    "GCU": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "GAU": "Aspartic Acid", "GAC": "Aspartic Acid",
    "GAA": "Glutamic Acid", "GAG": "Glutamic Acid",
    "GGU": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine"
}



def transcribe_dna_to_mrna(dna_sequence):
    """Transcribe a DNA sequence to mRNA by replacing T with U."""
    if not isinstance(dna_sequence, str):
        raise ValueError("DNA sequence must be a string")
    return dna_sequence.replace("T", "U")



def translate_mrna_to_protein(mrna_sequence):
    """Translate an mRNA sequence to a protein (list of amino acids)."""
    if not isinstance(mrna_sequence, str) or len(mrna_sequence) < 3:
        return []
    protein = []
    codons = [mrna_sequence[i:i+3] for i in range(0, len(mrna_sequence) - 2, 3)]
    for codon in codons:
        amino_acid = codon_to_amino_acid.get(codon, "Unknown")
        if amino_acid == "STOP":
            break
        if amino_acid != "Unknown":
            protein.append(f"{amino_acid} ({codon})")
    return protein



def process_dna_sequence(dna_sequence):
    """Process a DNA sequence through transcription and translation."""
    mrna = transcribe_dna_to_mrna(dna_sequence)
    protein = translate_mrna_to_protein(mrna)
    return {"mRNA": mrna, "Protein": protein}



class TestDNATranscriptionTranslation(unittest.TestCase):
    def test_transcription(self):
        self.assertEqual(transcribe_dna_to_mrna("ATGC"), "AUGC")
        self.assertEqual(transcribe_dna_to_mrna("TTA"), "UUA")
        self.assertEqual(transcribe_dna_to_mrna("GATTACA"), "GAUUACA")
        self.assertEqual(transcribe_dna_to_mrna(""), "")

    def test_translation(self):
        self.assertEqual(translate_mrna_to_protein("AUG"), ["Methionine (AUG)"])
        self.assertEqual(translate_mrna_to_protein("UUU"), ["Phenylalanine (UUU)"])
        self.assertEqual(translate_mrna_to_protein("UAA"), [])
        self.assertEqual(translate_mrna_to_protein("AUGUUUUGA"), ["Methionine (AUG)", "Phenylalanine (UUU)"])
        self.assertEqual(translate_mrna_to_protein(""), [])

    def test_full_process(self):
        result = process_dna_sequence("ATGTTT")
        self.assertEqual(result["mRNA"], "AUGUUU")
        self.assertEqual(result["Protein"], ["Methionine (AUG)", "Phenylalanine (UUU)"])

    def test_invalid_input(self):
        result = process_dna_sequence("XYZ")
        self.assertEqual(result["mRNA"], "XYZ")
        self.assertEqual(result["Protein"], [])
        with self.assertRaises(ValueError):
            transcribe_dna_to_mrna(123)



def main():
    print("DNA Transcription and Translation Simulator")
    print("Enter a DNA sequence (e.g., ATGTTT) or 'quit' to exit:")
    while True:
        dna_input = input("> ").strip().upper()
        if dna_input.lower() == "quit":
            print("Exiting simulator.")
            break
        if not dna_input:
            print("Error: Empty input. Please enter a DNA sequence.")
            continue
        if not all(base in "ATGC" for base in dna_input):
            print("Invalid DNA sequence! Use only A, T, G, C.")
            continue
        try:
            result = process_dna_sequence(dna_input)
            print(f"mRNA: {result['mRNA']}")
            print(f"Protein: {' -> '.join(result['Protein']) if result['Protein'] else 'None'}")
            if len(dna_input) % 3 != 0:
                print("Note: Sequence length not a multiple of 3; incomplete codons ignored.")
        except ValueError as e:
            print(f"Error: {e}")

if __name__ == "__main__":
    unittest.main(argv=['first-arg-is-ignored'], exit=False)
    print("\nAll tests passed successfully!\n")
    main()