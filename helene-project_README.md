```python
# Import necessary libraries
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Phylo, Entrez
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import AlignInfo
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from logomaker import Logo
from logomaker import alignment_to_matrix

```

## 1. Importing the Necessary Libraries

**From BioPython:**
- Bio.Blast: used to submit sequences to NCBI databases and search for homologues across species 
- Bio.SeqIO: used for parsing sequences from input files for alignment and writing alignments to output files in FASTA format
- Bio.Phylo: used for drawing Phylogenetic trees (Phylo.draw function) and (Bio.Phylo.TreeConstruction)
- Bio.Align: used for multiple and pairwise alignment of seqences, and (Align.Info) tool is used for extracting positon-specific scoring matricies (PSSM scores)
- Bio.SeqRecord: adds metadata such as descriptions and IDs to sequences
- Bio.Entrez: Connects to NCBI's Entrez database to retreieve sequence and domain data

**From matplotlib:**
- pyplot.figure: used to customize figure labels, size and layout

**From numpy:**
- np.array: used to covert data to a format which can be used in plots

**From collections:**
- Counter: used to summarize the amino acid composition across the alignment and count the frequency of amino acids in the alignmed sequences

**From Logomaker:**
- alignment_to_matrix: converts the alignment to a frequency matrix 
- Logo: used to convert the frequency matrix to generate a sequence logo




```python
# Performing a BLAST search to search for homoglogues
def perform_blast(query_file, output_file):
    print("Running BLAST...")
    with open(query_file) as f:
        blast_results = NCBIWWW.qblast("blastp", "nr", f.read())
    
    with open(output_file, "w") as out_file:
        out_file.write(blast_results.read())
    print(f"BLAST results saved to {output_file}")

```

## 2. Define Functions for Sequence Analysis

**Aim:** Search NCBI protein database for homogogous sequences to our queries

*perform_blast(query_file, output_file)*
- uses query in FASTA format
- performs BLAST search for homologues using (qblast)
- saves BLAST results in an XML file.



```python
# Parse BLAST results and save homologs in FASTA format
def parse_blast_results(blast_file, fasta_output, e_value_threshold=1e-5):
    print("Parsing BLAST results...")
    with open(blast_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        with open(fasta_output, "w") as output_handle:
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < e_value_threshold:
                            output_handle.write(f">{alignment.title}\n")
                            output_handle.write(f"{hsp.sbjct}\n")
    print(f"Homologous sequences saved to {fasta_output}")
```

**Aim:** to filter and save relavent homologues for further analysis 

*parse_blast_results(blast_file, fasta_output, e_value_threshold=1e-5)*
- Used to parse the BLAST XML results file
- Extrat the sequences with an e-value below the specified threshold (used 1e-5)
- Saves the homologous sequences to a FASTA file 


```python
# Retrieve additional homologous sequences
def retrieve_homologs(query_file, homologs_file):
    print("Retrieving homologous sequences...")
    Entrez.email = "helenerittscher@gmail.com"  # Replace with your email
    with open(query_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            query_id = record.id
            print(f"Searching homologs for: {query_id}")
            handle = Entrez.esearch(db="protein", term=query_id, retmax=10)
            record_ids = Entrez.read(handle)['IdList']
            with open(homologs_file, "w") as output_handle:
                for prot_id in record_ids:
                    prot_handle = Entrez.efetch(db="protein", id=prot_id, rettype="fasta", retmode="text")
                    output_handle.write(prot_handle.read())
    print(f"Homologous sequences saved to {homologs_file}")
```

**Aim:** To retrieve homologous protein sequences from Entrez. 
(This does not rely on BLAST so by querying this database it broadens the scope of the search and may reduce database biases as well as improve reproducability.)

*retrieve_homologs(query_file, homologs_file)*
- Uses Entrez API to seach for homologous sequences
- Retrieves homologous sequences in FASTA format and saves them 


```python
# Function to pad sequences to the same length
def pad_sequences(sequences):
    print("Padding sequences to the same length...")
    max_length = max(len(seq.seq) for seq in sequences)
    padded_sequences = []
    for seq in sequences:
        padded_seq = SeqRecord(Seq(str(seq.seq).ljust(max_length, "-")), id=seq.id, description=seq.description)
        padded_sequences.append(padded_seq)
    return padded_sequences
```

**Aim:** To ensure all sequences are aligned to the same length, so they can be used in the alignment tools in later steps

*pad_sequences(sequences)*
- Pads all sequences to the same length by appending gaps (-) where needed
- returns a list of padded sequences


```python
# Performing sequence alignment using Biopython PairwiseAligner
def align_sequences(input_file, output_file):
    print("Aligning sequences using Biopython PairwiseAligner...")
    aligner = PairwiseAligner()
    aligner.mode = "global"  # Use global alignment

    sequences = list(SeqIO.parse(input_file, "fasta"))
    if len(sequences) < 2:
        print("Error: At least two sequences are required for alignment.")
        return

    padded_sequences = pad_sequences(sequences)
    msa = MultipleSeqAlignment(padded_sequences)

    with open(output_file, "w") as output_handle:
        SeqIO.write(msa, output_handle, "fasta")

    print(f"Alignment saved to {output_file}")
```

**Aim:** To align sequences for phylogenetic and analysis of motifs

*align_sequences(input_file, output_file)*

- PairwiseAligner from Biopythin used for gloal alignmal on sequences from the FASTA input file
- pads sequences to the same length before alignment
- aligned sequences are saved to an output FASTA file



```python
# Building and visualizing a phylogenetic tree with better spacing
def build_phylogenetic_tree(alignment_file, tree_output):
    print("Building phylogenetic tree...")
    alignment = MultipleSeqAlignment(SeqIO.parse(alignment_file, "fasta"))

    # Calculate the distance matrix
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)

    # Construct the phylogenetic tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)

    # Save the tree in Newick format
    Phylo.write(tree, tree_output, "newick")
    print(f"Phylogenetic tree saved to {tree_output}")

    # Visualize the tree with advanced layout adjustments
    fig = plt.figure(figsize=(16, 12))  # Larger figure size for clarity
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)

    # Customizing the Phylogenetic Tree Figure
    plt.title("Phylogenetic Tree", fontsize=16)
    plt.xlabel("Branch Length", fontsize=12)
    plt.ylabel("Taxa", fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.tight_layout()
    plt.show()
```

**Aim:** Build a phylogenetic tree to analyze the evolutionary relationships between homologous sequences

*build_phylogenetic_tree(alignment_file, tree_output)*
- Calculates a distance matrix from the aigned sequences
- Uses the neighbor-joining mehod to contruct a phylogenetic tree
- Tree saved in Newick format and made into a customizable figure using Matplotlib


```python
# Function to analyze conserved domains using Entrez
# Requires valid Entrez email registration
def analyze_conserved_domains(fasta_file):
    print("Analyzing conserved domains...")
    Entrez.email = "helenerittscher@gmail.com"  # Replace with your email
    with open(fasta_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            print(f"Fetching domains for: {record.id}")
            handle = Entrez.efetch(db="protein", id=record.id, rettype="gp", retmode="text")
            print(handle.read())  # Modify to save or parse results

```

**Aim:** Identify conserved domains within the protein sequences

*analyze_conserved_domains(fasta_file)*
- Retrieves domain information for the protein sequences using Entrez API
- Prints domain data


```python
# Visualizing motifs 
def visualize_motifs(alignment_file):
    print("Visualizing motifs...")
    alignment = list(SeqIO.parse(alignment_file, "fasta"))

    # Convert alignment into a position-specific scoring matrix (PSSM)
    msa = MultipleSeqAlignment(alignment)
    summary_align = AlignInfo.SummaryInfo(msa)
    pssm = summary_align.pos_specific_score_matrix()

    # Truncate PSSM to the first 230 positions
    truncated_data = np.array([list(column.values())[:230] for column in pssm]).T

    # Visualize truncated PSSM as a heatmap
    plt.figure(figsize=(14, 8))  # Adjust figure size for clarity
    plt.imshow(truncated_data, aspect="auto", cmap="viridis")
    plt.colorbar(label="Score")
    plt.title("Motif Visualization (PSSM Heatmap)")
    plt.xlabel("Position in Alignment")
    plt.ylabel("Amino Acids")
    plt.xlim((0,230))
    plt.tight_layout()
    plt.show()
```

**Aim:** To identify and visualize the conserved motifs in the protein sequences

*visualize_motifs(alignment_file)*
- Uses the alignment to extract a position-specific scoring matrix (PSSM)
- A heatmap os the PSSM is created for the first 230 residues as sequences do not exceed this length
- This figure helps to visualize conserved motifs in the aligned sequences


```python
# Calculating amino acid frequencies
def compute_amino_acid_frequencies(alignment_file):
    print("Computing amino acid frequencies...")
    alignment = list(SeqIO.parse(alignment_file, "fasta"))
    concatenated_seqs = "".join([str(record.seq) for record in alignment])
    frequencies = Counter(concatenated_seqs)
    print("Amino acid frequencies:")
    for aa, freq in frequencies.items():
        print(f"{aa}: {freq}")
    return frequencies

```

**Aim:** Analyze the amino acd composion in the aligned sequences

*compute_amino_acid_frequencies(aligment_file)
- counts the frequency of each amino acid across all sequences in the alignment
- prints the frequencies for interpretation
  


```python
# Visualizing the sequence motifs using a sequence logo ( 230 residues only)
def visualize_sequence_logo(alignment_file):
    print("Visualizing sequence motifs using sequence logo...")
    alignment = list(SeqIO.parse(alignment_file, "fasta"))

    # Convert alignment into a frequency matrix
    sequences = [str(record.seq)[:230] for record in alignment]  # Truncate to first 200 residues
    alignment_matrix = alignment_to_matrix(sequences)

    # Plot the sequence logo
    plt.figure(figsize=(14, 8))  # Wider figure size for better spacing
    logo = Logo(alignment_matrix, font_name='Arial', color_scheme='chemistry')
    logo.style_spines(visible=False)
    logo.ax.set_xticks(range(0, 231, 20))  # Add ticks every 20 residues
    logo.ax.set_xticklabels(range(1, 232, 20), fontsize=10)  # Adjust tick font size
    logo.ax.set_title("Sequence Logo", fontsize=16)
    logo.ax.set_xlabel("Residue Position", fontsize=12)
    logo.ax.set_ylabel("Information Content", fontsize=12)
    plt.tight_layout()
    plt.show()
```

**Aim:** Visualize the conserved sequence features 

*Visualize_sequence_logo(alignment_file)*
- converts the sequence akignment into a frequency matrix for 230 residues
- uses Logomaker to create a sequence logo visualization
- Adjusts axis labels and formatting 


```python
# Main script
def main():
    # Input query sequences for BLAST
    query_file = "query_sequence.fasta"  # Input query sequence in FASTA format

    # Define sequences for RABC2A and RABD1
    sequences = {
        "RABC2A": "MGSSSGQSGYDLSFKILLIGDSGVGKSSLLVSFISSSVEDLAPTIGVDFKIKQLTVGGKRLKLTIWDTAGQERFRTLTSSYYRGAQGIILVYDVTRRETFTNLVDVWGKEIELYSTNQECVRMLVGNKVDRESERGVSREEGIALAKELNCMFLECSARTRQNVEQCFEELALKIMEVPSLLEEGSSAVKRNILKQKPEHQTNTQSGCCS",
        "RABD1": "MSNEYDYLFKLLLIGDSSVGKSCLLLRFADDAYIDSYISTIGVDFKIRTIEQDGKTIKLQIWDTAGQERFRTITSSYYRGAHGIIIVYDCTEMESFNNVKQWLSEIDRYANESVCKLLIGNKNDMVESKVVSTETGRALADELGIPFLETSAKDSINVEQAFLTIAGEIKKKMGSQTNANKTSGPGTVQMKGQPIQQNNGGCCGQ"
    }

    # Write sequences to the query file
    with open(query_file, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")

    homologs_fasta = "homologs.fasta"  # Homologous sequences
    blast_output = "blast_results.xml"  # BLAST results (optional if using BLAST)
    aligned_file = "aligned_homologs.fasta"  # Aligned sequences
    tree_file = "aligned_homologs.dnd"  # Phylogenetic tree

    # Step 1: Retrieve homologous sequences
    retrieve_homologs(query_file, homologs_fasta)

    # Step 2: Align sequences
    align_sequences(homologs_fasta, aligned_file)

    # Step 3: Build and visualize a phylogenetic tree
    build_phylogenetic_tree(aligned_file, tree_file)

    # Step 4: Analyze conserved domains
    analyze_conserved_domains(homologs_fasta)

    # Step 5: Visualize motifs as heatmap
    visualize_motifs(aligned_file)

    # Step 6: Visualize motifs as sequence logo
    visualize_sequence_logo(aligned_file)

if __name__ == "__main__":
    main()
```

    Retrieving homologous sequences...
    Searching homologs for: RABC2A
    Searching homologs for: RABD1
    Homologous sequences saved to homologs.fasta
    Aligning sequences using Biopython PairwiseAligner...
    Padding sequences to the same length...
    Alignment saved to aligned_homologs.fasta
    Building phylogenetic tree...
    Phylogenetic tree saved to aligned_homologs.dnd



    
![png](output_22_1.png)
    


    Analyzing conserved domains...
    Fetching domains for: KAL3651522.1
    LOCUS       KAL3651522               203 aa            linear   PLN 17-DEC-2024
    DEFINITION  Ras-related protein rabd1 [Castilleja foliolosa].
    ACCESSION   KAL3651522
    VERSION     KAL3651522.1
    DBLINK      BioProject: PRJNA1011274
                BioSample: SAMN37213163
    DBSOURCE    accession JAVIJP010000006.1
    KEYWORDS    .
    SOURCE      Castilleja foliolosa
      ORGANISM  Castilleja foliolosa
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Lamiales; Orobanchaceae;
                Pedicularideae; Castillejinae; Castilleja.
    REFERENCE   1  (residues 1 to 203)
      AUTHORS   Buerger,M., Peterson,D. and Chory,J.
      TITLE     Strigolactones Initiate the Formation of Haustorium-like Structures
                in Castilleja
      JOURNAL   iScience, 111491 (2024)
      REMARK    DOI: 10.1016/j.isci.2024.111491
    REFERENCE   2  (residues 1 to 203)
      AUTHORS   Burger,M. and Chory,J.
      TITLE     Direct Submission
      JOURNAL   Submitted (23-NOV-2024) Plant Biology Department, Salk Institute
                for Biological Studies, 10010 N Torrey Pines Rd, La Jolla, CA
                92037, USA
    COMMENT     ##Genome-Assembly-Data-START##
                Assembly Method        :: Hifiasm v. 0.19.1
                Genome Representation  :: Full
                Expected Final Version :: Yes
                Genome Coverage        :: 15.0x
                Sequencing Technology  :: PacBio Sequel
                ##Genome-Assembly-Data-END##
                Method: conceptual translation.
    FEATURES             Location/Qualifiers
         source          1..203
                         /organism="Castilleja foliolosa"
                         /isolate="Tecolote"
                         /db_xref="taxon:1961234"
                         /tissue_type="flower"
                         /dev_stage="Adult plant"
                         /geo_loc_name="USA: San Diego County, CA"
                         /collection_date="2022-01-31"
         Protein         1..203
                         /product="Ras-related protein rabd1"
         CDS             1..203
                         /gene="RABD1_1"
                         /locus_tag="CASFOL_004524"
                         /coded_by="join(JAVIJP010000006.1:3572658..3572671,
                         JAVIJP010000006.1:3572759..3572831,
                         JAVIJP010000006.1:3573112..3573159,
                         JAVIJP010000006.1:3573265..3573312,
                         JAVIJP010000006.1:3573656..3573727,
                         JAVIJP010000006.1:3574006..3574161,
                         JAVIJP010000006.1:3574291..3574394,
                         JAVIJP010000006.1:3574651..3574747)"
                         /note="COG:U; EggNog:ENOG503GD0A;
                         GO_function: GO:0003924 - GTPase activity [Evidence IEA];
                         GO_function: GO:0005525 - GTP binding [Evidence IEA]"
                         /db_xref="InterPro:IPR001806"
                         /db_xref="InterPro:IPR005225"
                         /db_xref="InterPro:IPR027417"
                         /db_xref="PFAM:PF00025"
                         /db_xref="PFAM:PF00071"
                         /db_xref="PFAM:PF01926"
                         /db_xref="PFAM:PF04670"
                         /db_xref="PFAM:PF08477"
    ORIGIN      
            1 mgseydylfk llligdssvg ksclllrfad dsyvdsyist igvdfkirtv eldgktiklq
           61 iwdtagqerf rtitssyyrg ahgiiivydv temesfnnvk qwlseidrya sssvckllvg
          121 nkcdlvdskv vdtqtgkala delgipflet sakdsinveq afltmageik kksgnqpman
          181 kaptgrvqis gqpieqnssn ccg
    //
    
    
    Fetching domains for: KAL3627834.1
    LOCUS       KAL3627834               149 aa            linear   PLN 17-DEC-2024
    DEFINITION  Ras-related protein rabd1 [Castilleja foliolosa].
    ACCESSION   KAL3627834
    VERSION     KAL3627834.1
    DBLINK      BioProject: PRJNA1011274
                BioSample: SAMN37213163
    DBSOURCE    accession JAVIJP010000038.1
    KEYWORDS    .
    SOURCE      Castilleja foliolosa
      ORGANISM  Castilleja foliolosa
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Lamiales; Orobanchaceae;
                Pedicularideae; Castillejinae; Castilleja.
    REFERENCE   1  (residues 1 to 149)
      AUTHORS   Buerger,M., Peterson,D. and Chory,J.
      TITLE     Strigolactones Initiate the Formation of Haustorium-like Structures
                in Castilleja
      JOURNAL   iScience, 111491 (2024)
      REMARK    DOI: 10.1016/j.isci.2024.111491
    REFERENCE   2  (residues 1 to 149)
      AUTHORS   Burger,M. and Chory,J.
      TITLE     Direct Submission
      JOURNAL   Submitted (23-NOV-2024) Plant Biology Department, Salk Institute
                for Biological Studies, 10010 N Torrey Pines Rd, La Jolla, CA
                92037, USA
    COMMENT     ##Genome-Assembly-Data-START##
                Assembly Method        :: Hifiasm v. 0.19.1
                Genome Representation  :: Full
                Expected Final Version :: Yes
                Genome Coverage        :: 15.0x
                Sequencing Technology  :: PacBio Sequel
                ##Genome-Assembly-Data-END##
                Method: conceptual translation.
    FEATURES             Location/Qualifiers
         source          1..149
                         /organism="Castilleja foliolosa"
                         /isolate="Tecolote"
                         /db_xref="taxon:1961234"
                         /tissue_type="flower"
                         /dev_stage="Adult plant"
                         /geo_loc_name="USA: San Diego County, CA"
                         /collection_date="2022-01-31"
         Protein         1..149
                         /product="Ras-related protein rabd1"
         CDS             1..149
                         /gene="RABD1_2"
                         /locus_tag="CASFOL_028249"
                         /coded_by="join(JAVIJP010000038.1:2236425..2236513,
                         JAVIJP010000038.1:2236727..2236733,
                         JAVIJP010000038.1:2236940..2237095,
                         JAVIJP010000038.1:2237325..2237428,
                         JAVIJP010000038.1:2237839..2237927,
                         JAVIJP010000038.1:2239759..2239763)"
                         /note="COG:U; EggNog:ENOG503GD0A;
                         GO_function: GO:0003924 - GTPase activity [Evidence IEA];
                         GO_function: GO:0005525 - GTP binding [Evidence IEA]"
                         /db_xref="InterPro:IPR001806"
                         /db_xref="InterPro:IPR005225"
                         /db_xref="InterPro:IPR027417"
                         /db_xref="PFAM:PF00071"
                         /db_xref="PFAM:PF08477"
    ORIGIN      
            1 mgfkwdtagq erfrtitssy yrgahgmasf nmivydvtem esfnnvkqwl neidryands
           61 vckllvgnkc dlveskvvdt qtakefadel gipfletsak dainveqafl tmageikkkm
          121 gnqptgnrks gstvqikgqp ieqksnccg
    //
    
    
    Fetching domains for: KAL3620804.1
    LOCUS       KAL3620804               202 aa            linear   PLN 17-DEC-2024
    DEFINITION  Ras-related protein rabd1 [Castilleja foliolosa].
    ACCESSION   KAL3620804
    VERSION     KAL3620804.1
    DBLINK      BioProject: PRJNA1011274
                BioSample: SAMN37213163
    DBSOURCE    accession JAVIJP010000066.1
    KEYWORDS    .
    SOURCE      Castilleja foliolosa
      ORGANISM  Castilleja foliolosa
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Lamiales; Orobanchaceae;
                Pedicularideae; Castillejinae; Castilleja.
    REFERENCE   1  (residues 1 to 202)
      AUTHORS   Buerger,M., Peterson,D. and Chory,J.
      TITLE     Strigolactones Initiate the Formation of Haustorium-like Structures
                in Castilleja
      JOURNAL   iScience, 111491 (2024)
      REMARK    DOI: 10.1016/j.isci.2024.111491
    REFERENCE   2  (residues 1 to 202)
      AUTHORS   Burger,M. and Chory,J.
      TITLE     Direct Submission
      JOURNAL   Submitted (23-NOV-2024) Plant Biology Department, Salk Institute
                for Biological Studies, 10010 N Torrey Pines Rd, La Jolla, CA
                92037, USA
    COMMENT     ##Genome-Assembly-Data-START##
                Assembly Method        :: Hifiasm v. 0.19.1
                Genome Representation  :: Full
                Expected Final Version :: Yes
                Genome Coverage        :: 15.0x
                Sequencing Technology  :: PacBio Sequel
                ##Genome-Assembly-Data-END##
                Method: conceptual translation.
    FEATURES             Location/Qualifiers
         source          1..202
                         /organism="Castilleja foliolosa"
                         /isolate="Tecolote"
                         /db_xref="taxon:1961234"
                         /tissue_type="flower"
                         /dev_stage="Adult plant"
                         /geo_loc_name="USA: San Diego County, CA"
                         /collection_date="2022-01-31"
         Protein         1..202
                         /product="Ras-related protein rabd1"
         CDS             1..202
                         /gene="RABD1_3"
                         /locus_tag="CASFOL_035716"
                         /coded_by="join(JAVIJP010000066.1:17252629..17252642,
                         JAVIJP010000066.1:17252743..17252815,
                         JAVIJP010000066.1:17252952..17252999,
                         JAVIJP010000066.1:17253096..17253143,
                         JAVIJP010000066.1:17253715..17253786,
                         JAVIJP010000066.1:17254768..17254923,
                         JAVIJP010000066.1:17255125..17255228,
                         JAVIJP010000066.1:17255516..17255609)"
                         /note="COG:U; EggNog:ENOG503GD0A;
                         GO_function: GO:0003924 - GTPase activity [Evidence IEA];
                         GO_function: GO:0005525 - GTP binding [Evidence IEA]"
                         /db_xref="InterPro:IPR001806"
                         /db_xref="InterPro:IPR005225"
                         /db_xref="InterPro:IPR027417"
                         /db_xref="PFAM:PF00009"
                         /db_xref="PFAM:PF00025"
                         /db_xref="PFAM:PF00071"
                         /db_xref="PFAM:PF01926"
                         /db_xref="PFAM:PF04670"
                         /db_xref="PFAM:PF08477"
    ORIGIN      
            1 msseydylfk llligdssvg ksclllrfad dsyvdsyist igvdfkirtv eldgktsklq
           61 iwdtagqerf rtitssyyrg ahgiiivydv temesfnnvk qwlneidrya ndsvckllvg
          121 nkcdlveskv vdtqtakafa delgipflet sakdainveq afltmageik kkmsnqptgs
          181 kksastvqik gqpieqksnc cg
    //
    
    
    Fetching domains for: KAL3614467.1
    LOCUS       KAL3614467               110 aa            linear   PLN 17-DEC-2024
    DEFINITION  Ras-related protein rabd1 [Castilleja foliolosa].
    ACCESSION   KAL3614467
    VERSION     KAL3614467.1
    DBLINK      BioProject: PRJNA1011274
                BioSample: SAMN37213163
    DBSOURCE    accession JAVIJP010000105.1
    KEYWORDS    .
    SOURCE      Castilleja foliolosa
      ORGANISM  Castilleja foliolosa
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Lamiales; Orobanchaceae;
                Pedicularideae; Castillejinae; Castilleja.
    REFERENCE   1  (residues 1 to 110)
      AUTHORS   Buerger,M., Peterson,D. and Chory,J.
      TITLE     Strigolactones Initiate the Formation of Haustorium-like Structures
                in Castilleja
      JOURNAL   iScience, 111491 (2024)
      REMARK    DOI: 10.1016/j.isci.2024.111491
    REFERENCE   2  (residues 1 to 110)
      AUTHORS   Burger,M. and Chory,J.
      TITLE     Direct Submission
      JOURNAL   Submitted (23-NOV-2024) Plant Biology Department, Salk Institute
                for Biological Studies, 10010 N Torrey Pines Rd, La Jolla, CA
                92037, USA
    COMMENT     ##Genome-Assembly-Data-START##
                Assembly Method        :: Hifiasm v. 0.19.1
                Genome Representation  :: Full
                Expected Final Version :: Yes
                Genome Coverage        :: 15.0x
                Sequencing Technology  :: PacBio Sequel
                ##Genome-Assembly-Data-END##
                Method: conceptual translation.
    FEATURES             Location/Qualifiers
         source          1..110
                         /organism="Castilleja foliolosa"
                         /isolate="Tecolote"
                         /db_xref="taxon:1961234"
                         /tissue_type="flower"
                         /dev_stage="Adult plant"
                         /geo_loc_name="USA: San Diego County, CA"
                         /collection_date="2022-01-31"
         Protein         1..110
                         /product="Ras-related protein rabd1"
         CDS             1..110
                         /gene="RABD1_4"
                         /locus_tag="CASFOL_041553"
                         /coded_by="join(JAVIJP010000105.1:852292..852426,
                         JAVIJP010000105.1:852655..852758,
                         JAVIJP010000105.1:853169..853257,
                         JAVIJP010000105.1:855096..855100)"
                         /note="COG:U; EggNog:ENOG503GD0A;
                         GO_function: GO:0003924 - GTPase activity [Evidence IEA];
                         GO_function: GO:0005525 - GTP binding [Evidence IEA]"
                         /db_xref="InterPro:IPR001806"
                         /db_xref="InterPro:IPR027417"
    ORIGIN      
            1 mesfnnvkqw lneidryand svckllvgnk cdlveskvvd tqtakefade lgipfletsa
           61 kdainveqaf ltmageikkk mgnqptgnrk sgstvqikgq pieqksnccg
    //
    
    
    Fetching domains for: sp|Q9ZRE2.1|RABD1_ARATH
    LOCUS       RABD1_ARATH              205 aa            linear   PLN 27-NOV-2024
    DEFINITION  RecName: Full=Ras-related protein RABD1; Short=AtRABD1; AltName:
                Full=Ras-related protein ATFP8.
    ACCESSION   Q9ZRE2
    VERSION     Q9ZRE2.1
    DBSOURCE    UniProtKB: locus RABD1_ARATH, accession Q9ZRE2;
                class: standard.
                created: Sep 2, 2008.
                sequence updated: May 1, 1999.
                annotation updated: Nov 27, 2024.
                xrefs: U64911.1, AAD00111.1, AC016795.6, AAF23189.1, CP002686.1,
                AEE75090.1, AK118150.1, BAC42775.1, BT005576.1, AAO63996.1,
                NP_187779.1
                xrefs (non-sequence databases): AlphaFoldDB:Q9ZRE2, SMR:Q9ZRE2,
                BioGRID:5679, IntAct:Q9ZRE2, STRING:3702.Q9ZRE2, iPTMnet:Q9ZRE2,
                SwissPalm:Q9ZRE2, PaxDb:3702-AT3G11730.1, ProteomicsDB:236556,
                EnsemblPlants:AT3G11730.1, EnsemblPlants:AT3G11730.1,
                EnsemblPlants:AT3G11730, GeneID:820345, Gramene:AT3G11730.1,
                KEGG:ath:AT3G11730, Araport:AT3G11730, TAIR:AT3G11730,
                eggNOG:KOG0084, HOGENOM:CLU_041217_10_1_1, InParanoid:Q9ZRE2,
                OMA:QQNSNCC, OrthoDB:8685at2759, PhylomeDB:Q9ZRE2, PRO:PR:Q9ZRE2,
                Proteomes:UP000006548, ExpressionAtlas:Q9ZRE2, GO:0000139,
                GO:0005886, GO:0032588, GO:0005525, GO:0030742, GO:0003924,
                GO:0080115, GO:0006888, GO:0015031, CDD:cd01869,
                FunFam:3.40.50.300:FF:000069, Gene3D:3.40.50.300,
                InterPro:IPR027417, InterPro:IPR050227, InterPro:IPR005225,
                InterPro:IPR001806, NCBIfam:TIGR00231, PANTHER:PTHR47977,
                PANTHER:PTHR47977:SF22, Pfam:PF00071, PRINTS:PR00449,
                SMART:SM00177, SMART:SM00175, SMART:SM00176, SMART:SM00173,
                SMART:SM00174, SUPFAM:SSF52540, PROSITE:PS51419
    KEYWORDS    Acetylation; ER-Golgi transport; Golgi apparatus; GTP-binding;
                Lipoprotein; Membrane; Nucleotide-binding; Prenylation; Protein
                transport; Reference proteome; Transport.
    SOURCE      Arabidopsis thaliana (thale cress)
      ORGANISM  Arabidopsis thaliana
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; rosids; malvids; Brassicales; Brassicaceae;
                Camelineae; Arabidopsis.
    REFERENCE   1  (residues 1 to 205)
      AUTHORS   Biermann,B., Randall,S.K. and Crowell,D.N.
      TITLE     Identification and isoprenylation of plant GTP-binding proteins
      JOURNAL   Plant Mol Biol 31 (5), 1021-1028 (1996)
       PUBMED   8843944
      REMARK    NUCLEOTIDE SEQUENCE [MRNA], AND ISOPRENYLATION.
    REFERENCE   2  (residues 1 to 205)
      AUTHORS   Salanoubat,M., Lemcke,K., Rieger,M., Ansorge,W., Unseld,M.,
                Fartmann,B., Valle,G., Blocker,H., Perez-Alonso,M., Obermaier,B.,
                Delseny,M., Boutry,M., Grivell,L.A., Mache,R., Puigdomenech,P., De
                Simone,V., Choisne,N., Artiguenave,F., Robert,C., Brottier,P.,
                Wincker,P., Cattolico,L., Weissenbach,J., Saurin,W., Quetier,F.,
                Schafer,M., Muller-Auer,S., Gabel,C., Fuchs,M., Benes,V.,
                Wurmbach,E., Drzonek,H., Erfle,H., Jordan,N., Bangert,S.,
                Wiedelmann,R., Kranz,H., Voss,H., Holland,R., Brandt,P.,
                Nyakatura,G., Vezzi,A., D'Angelo,M., Pallavicini,A., Toppo,S.,
                Simionati,B., Conrad,A., Hornischer,K., Kauer,G., Lohnert,T.H.,
                Nordsiek,G., Reichelt,J., Scharfe,M., Schon,O., Bargues,M.,
                Terol,J., Climent,J., Navarro,P., Collado,C., Perez-Perez,A.,
                Ottenwalder,B., Duchemin,D., Cooke,R., Laudie,M., Berger-Llauro,C.,
                Purnelle,B., Masuy,D., de Haan,M., Maarse,A.C., Alcaraz,J.P.,
                Cottet,A., Casacuberta,E., Monfort,A., Argiriou,A., flores,M.,
                Liguori,R., Vitale,D., Mannhaupt,G., Haase,D., Schoof,H., Rudd,S.,
                Zaccaria,P., Mewes,H.W., Mayer,K.F., Kaul,S., Town,C.D., Koo,H.L.,
                Tallon,L.J., Jenkins,J., Rooney,T., Rizzo,M., Walts,A.,
                Utterback,T., Fujii,C.Y., Shea,T.P., Creasy,T.H., Haas,B.,
                Maiti,R., Wu,D., Peterson,J., Van Aken,S., Pai,G., Militscher,J.,
                Sellers,P., Gill,J.E., Feldblyum,T.V., Preuss,D., Lin,X.,
                Nierman,W.C., Salzberg,S.L., White,O., Venter,J.C., Fraser,C.M.,
                Kaneko,T., Nakamura,Y., Sato,S., Kato,T., Asamizu,E., Sasamoto,S.,
                Kimura,T., Idesawa,K., Kawashima,K., Kishida,Y., Kiyokawa,C.,
                Kohara,M., Matsumoto,M., Matsuno,A., Muraki,A., Nakayama,S.,
                Nakazaki,N., Shinpo,S., Takeuchi,C., Wada,T., Watanabe,A.,
                Yamada,M., Yasuda,M. and Tabata,S.
      CONSRTM   European Union Chromosome 3 Arabidopsis Sequencing Consortium;
                Institute for Genomic Research; Kazusa DNA Research Institute
      TITLE     Sequence and analysis of chromosome 3 of the plant Arabidopsis
                thaliana
      JOURNAL   Nature 408 (6814), 820-822 (2000)
       PUBMED   11130713
      REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].;
                STRAIN=cv. Columbia
    REFERENCE   3  (residues 1 to 205)
      AUTHORS   Cheng,C.Y., Krishnakumar,V., Chan,A.P., Thibaud-Nissen,F.,
                Schobel,S. and Town,C.D.
      TITLE     Araport11: a complete reannotation of the Arabidopsis thaliana
                reference genome
      JOURNAL   Plant J 89 (4), 789-804 (2017)
       PUBMED   27862469
      REMARK    GENOME REANNOTATION.;
                STRAIN=cv. Columbia
    REFERENCE   4  (residues 1 to 205)
      AUTHORS   Seki,M., Narusaka,M., Kamiya,A., Ishida,J., Satou,M., Sakurai,T.,
                Nakajima,M., Enju,A., Akiyama,K., Oono,Y., Muramatsu,M.,
                Hayashizaki,Y., Kawai,J., Carninci,P., Itoh,M., Ishii,Y.,
                Arakawa,T., Shibata,K., Shinagawa,A. and Shinozaki,K.
      TITLE     Functional annotation of a full-length Arabidopsis cDNA collection
      JOURNAL   Science 296 (5565), 141-145 (2002)
       PUBMED   11910074
      REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA].;
                STRAIN=cv. Columbia
    REFERENCE   5  (residues 1 to 205)
      AUTHORS   Yamada,K., Lim,J., Dale,J.M., Chen,H., Shinn,P., Palm,C.J.,
                Southwick,A.M., Wu,H.C., Kim,C., Nguyen,M., Pham,P., Cheuk,R.,
                Karlin-Newmann,G., Liu,S.X., Lam,B., Sakano,H., Wu,T., Yu,G.,
                Miranda,M., Quach,H.L., Tripp,M., Chang,C.H., Lee,J.M., Toriumi,M.,
                Chan,M.M., Tang,C.C., Onodera,C.S., Deng,J.M., Akiyama,K.,
                Ansari,Y., Arakawa,T., Banh,J., Banno,F., Bowser,L., Brooks,S.,
                Carninci,P., Chao,Q., Choy,N., Enju,A., Goldsmith,A.D., Gurjal,M.,
                Hansen,N.F., Hayashizaki,Y., Johnson-Hopson,C., Hsuan,V.W.,
                Iida,K., Karnes,M., Khan,S., Koesema,E., Ishida,J., Jiang,P.X.,
                Jones,T., Kawai,J., Kamiya,A., Meyers,C., Nakajima,M., Narusaka,M.,
                Seki,M., Sakurai,T., Satou,M., Tamse,R., Vaysberg,M.,
                Wallender,E.K., Wong,C., Yamamura,Y., Yuan,S., Shinozaki,K.,
                Davis,R.W., Theologis,A. and Ecker,J.R.
      TITLE     Empirical analysis of transcriptional activity in the Arabidopsis
                genome
      JOURNAL   Science 302 (5646), 842-846 (2003)
       PUBMED   14593172
      REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA].;
                STRAIN=cv. Columbia
    REFERENCE   6  (residues 1 to 205)
      AUTHORS   Vernoud,V., Horton,A.C., Yang,Z. and Nielsen,E.
      TITLE     Analysis of the small GTPase gene superfamily of Arabidopsis
      JOURNAL   Plant Physiol 131 (3), 1191-1208 (2003)
       PUBMED   12644670
      REMARK    GENE FAMILY, AND NOMENCLATURE.
    REFERENCE   7  (residues 1 to 205)
      AUTHORS   Latijnhouwers,M., Gillespie,T., Boevink,P., Kriechbaumer,V.,
                Hawes,C. and Carvalho,C.M.
      TITLE     Localization and domain characterization of Arabidopsis golgin
                candidates
      JOURNAL   J Exp Bot 58 (15-16), 4373-4386 (2007)
       PUBMED   18182439
      REMARK    INTERACTION WITH GC5.
    REFERENCE   8  (residues 1 to 205)
      AUTHORS   Hashimoto,K., Igarashi,H., Mano,S., Takenaka,C., Shiina,T.,
                Yamaguchi,M., Demura,T., Nishimura,M., Shimmen,T. and Yokota,E.
      TITLE     An isoform of Arabidopsis myosin XI interacts with small GTPases in
                its C-terminal tail region
      JOURNAL   J Exp Bot 59 (13), 3523-3531 (2008)
       PUBMED   18703495
      REMARK    INTERACTION WITH XI-2/MYA2.
    REFERENCE   9  (residues 1 to 205)
      AUTHORS   Pinheiro,H., Samalova,M., Geldner,N., Chory,J., Martinez,A. and
                Moore,I.
      TITLE     Genetic evidence that the higher plant Rab-D1 and Rab-D2 GTPases
                exhibit distinct but overlapping interactions in the early
                secretory pathway
      JOURNAL   J Cell Sci 122 (Pt 20), 3749-3758 (2009)
       PUBMED   19789181
      REMARK    SUBCELLULAR LOCATION, AND MUTAGENESIS OF SER-22; GLN-67 AND
                ASN-121.
    REFERENCE   10 (residues 1 to 205)
      AUTHORS   Bienvenut,W.V., Sumpton,D., Martinez,A., Lilla,S., Espagne,C.,
                Meinnel,T. and Giglione,C.
      TITLE     Comparative large scale characterization of plant versus mammal
                proteins reveals similar and idiosyncratic N-alpha-acetylation
                features
      JOURNAL   Mol Cell Proteomics 11 (6), M111.015131 (2012)
       PUBMED   22223895
      REMARK    ACETYLATION [LARGE SCALE ANALYSIS] AT SER-2, CLEAVAGE OF INITIATOR
                METHIONINE [LARGE SCALE ANALYSIS], AND IDENTIFICATION BY MASS
                SPECTROMETRY [LARGE SCALE ANALYSIS].
    COMMENT     [FUNCTION] Protein transport. Regulator of membrane traffic from
                the Golgi apparatus towards the endoplasmic reticulum (ER).
                [SUBUNIT] Does not interact with GC5. Interacts with XI-2/MYA2.
                {ECO:0000269|PubMed:18182439, ECO:0000269|PubMed:18703495}.
                [INTERACTION] Q9ZRE2; Q9LKB9: XI-2; NbExp=3; IntAct=EBI-2009542,
                EBI-2009528.
                [SUBCELLULAR LOCATION] Golgi apparatus, trans-Golgi network
                membrane {ECO:0000269|PubMed:19789181}. Golgi apparatus membrane
                {ECO:0000305|PubMed:19789181}; Lipid-anchor
                {ECO:0000305|PubMed:19789181}.
                [SIMILARITY] Belongs to the small GTPase superfamily. Rab family.
                {ECO:0000305}.
    FEATURES             Location/Qualifiers
         source          1..205
                         /organism="Arabidopsis thaliana"
                         /db_xref="taxon:3702"
         gene            1..205
                         /gene="RABD1"
                         /locus_tag="At3g11730"
                         /gene_synonym="ATFP8"
                         /gene_synonym="F26K24.2"
         Protein         1..205
                         /product="Ras-related protein RABD1"
                         /note="AtRABD1; Ras-related protein ATFP8"
                         /UniProtKB_evidence="Evidence at protein level"
         Region          2..205
                         /region_name="Mature chain"
                         /note="Ras-related protein RABD1. /id=PRO_0000348544."
         Site            2
                         /site_type="acetylation"
                         /note="N-acetylserine.
                         /evidence=ECO:0007744|PubMed:22223895."
         Region          7..172
                         /region_name="Rab1_Ypt1"
                         /note="Rab GTPase family 1 includes the yeast homolog
                         Ypt1; cd01869"
                         /db_xref="CDD:206661"
         Site            7..10
                         /site_type="other"
                         /note="Rab subfamily motif 1 (RabSF1)"
                         /db_xref="CDD:206661"
         Site            15..22
                         /site_type="other"
                         /note="G1 box"
                         /db_xref="CDD:206661"
         Site            order(16,18..23,33,121..122,124,151..153)
                         /site_type="other"
                         /note="GTP/Mg2+ binding site [chemical binding]"
                         /db_xref="CDD:206661"
         Site            22
                         /site_type="mutagenized"
                         /note="S->N: Dominant negative (GDP-bound form); no effect
                         on trafficking between the ER and Golgi.
                         /evidence=ECO:0000269|PubMed:19789181."
         Site            order(23..34,37..39)
                         /site_type="other"
                         /note="Rab subfamily motif 2 (RabSF2)"
                         /db_xref="CDD:206661"
         Site            order(33,38..46)
                         /site_type="other"
                         /note="Switch I region"
                         /db_xref="CDD:206661"
         Region          37..45
                         /region_name="Short sequence motif of biological interest"
                         /note="Effector region. /evidence=ECO:0000250."
         Site            order(38,42..49,56,58)
                         /site_type="other"
                         /note="putative GEF interaction site [polypeptide
                         binding]"
                         /db_xref="CDD:206661"
         Site            40
                         /site_type="other"
                         /note="G2 box"
                         /db_xref="CDD:206661"
         Site            order(41,43..45,60,62,69..70,73,77,79..82,169)
                         /site_type="active"
                         /note="putative effector interaction site [active]"
                         /db_xref="CDD:206661"
         Site            order(41..42,44,62..63,70,72,74..76,79)
                         /site_type="active"
                         /note="GDI interaction site [active]"
                         /db_xref="CDD:206661"
         Site            41..45
                         /site_type="other"
                         /note="Rab family motif 1 (RabF1)"
                         /db_xref="CDD:206661"
         Site            58..62
                         /site_type="other"
                         /note="Rab family motif 2 (RabF2)"
                         /db_xref="CDD:206661"
         Site            63..66
                         /site_type="other"
                         /note="G3 box"
                         /db_xref="CDD:206661"
         Site            order(66,68..78)
                         /site_type="other"
                         /note="Switch II region"
                         /db_xref="CDD:206661"
         Site            67
                         /site_type="mutagenized"
                         /note="Q->L: No effect on trafficking between the ER and
                         Golgi. /evidence=ECO:0000269|PubMed:19789181."
         Site            69..74
                         /site_type="other"
                         /note="Rab family motif 3 (RabF3)"
                         /db_xref="CDD:206661"
         Site            77..81
                         /site_type="other"
                         /note="Rab family motif 4 (RabF4)"
                         /db_xref="CDD:206661"
         Site            86..91
                         /site_type="other"
                         /note="Rab family motif 5 (RabF5)"
                         /db_xref="CDD:206661"
         Site            112..118
                         /site_type="other"
                         /note="Rab subfamily motif 3 (RabSF3)"
                         /db_xref="CDD:206661"
         Site            121..124
                         /site_type="other"
                         /note="G4 box"
                         /db_xref="CDD:206661"
         Site            121
                         /site_type="mutagenized"
                         /note="N->I: Inhibits trafficking between the ER and
                         Golgi. Causes severe dwarfism and seedlings death.
                         /evidence=ECO:0000269|PubMed:19789181."
         Site            151..153
                         /site_type="other"
                         /note="G5 box"
                         /db_xref="CDD:206661"
         Site            167..172
                         /site_type="other"
                         /note="Rab subfamily motif 4 (RabSF4)"
                         /db_xref="CDD:206661"
         Region          174..205
                         /region_name="Region of interest in the sequence"
                         /note="Disordered. /evidence=ECO:0000256|SAM:MobiDB-lite."
         Site            202
                         /site_type="lipid-binding"
                         /note="S-geranylgeranyl cysteine. /evidence=ECO:0000250."
         Site            203
                         /site_type="lipid-binding"
                         /note="S-geranylgeranyl cysteine. /evidence=ECO:0000250."
    ORIGIN      
            1 msneydylfk llligdssvg ksclllrfad dayidsyist igvdfkirti eqdgktiklq
           61 iwdtagqerf rtitssyyrg ahgiiivydc temesfnnvk qwlseidrya nesvckllig
          121 nkndmveskv vstetgrala delgipflet sakdsinveq afltiageik kkmgsqtnan
          181 ktsgpgtvqm kgqpiqqnng gccgq
    //
    
    
    Fetching domains for: sp|Q9LKB9.1|MYO6_ARATH
    LOCUS       MYO6_ARATH              1505 aa            linear   PLN 27-NOV-2024
    DEFINITION  RecName: Full=Myosin-6; AltName: Full=AtMYA2.
    ACCESSION   Q9LKB9
    VERSION     Q9LKB9.1
    DBSOURCE    UniProtKB: locus MYO6_ARATH, accession Q9LKB9;
                class: standard.
                extra accessions:Q0WPF0,Q39158
                created: Jun 26, 2013.
                sequence updated: Oct 1, 2000.
                annotation updated: Nov 27, 2024.
                xrefs: Z34293.1, CAA84066.1, AP000368.1, BAA98070.1, CP002688.1,
                AED95022.1, AK229124.1, BAF00999.1, S51824, NP_199203.1, 7DHW_A
                xrefs (non-sequence databases): PDBsum:7DHW, AlphaFoldDB:Q9LKB9,
                SMR:Q9LKB9, BioGRID:19662, IntAct:Q9LKB9, STRING:3702.Q9LKB9,
                iPTMnet:Q9LKB9, PaxDb:3702-AT5G43900.3, ProteomicsDB:251405,
                EnsemblPlants:AT5G43900.1, EnsemblPlants:AT5G43900.1,
                EnsemblPlants:AT5G43900, GeneID:834412, Gramene:AT5G43900.1,
                KEGG:ath:AT5G43900, Araport:AT5G43900, TAIR:AT5G43900,
                eggNOG:KOG0160, HOGENOM:CLU_000192_3_1_1, InParanoid:Q9LKB9,
                OMA:NQDNNDH, PhylomeDB:Q9LKB9, PRO:PR:Q9LKB9,
                Proteomes:UP000006548, ExpressionAtlas:Q9LKB9, GO:0015629,
                GO:0005737, GO:0016020, GO:0016459, GO:0031982, GO:0051015,
                GO:0005524, GO:0005516, GO:0000146, GO:0007015, GO:0048767,
                GO:0030050, CDD:cd15475, CDD:cd01384, FunFam:1.20.58.530:FF:000002,
                FunFam:1.20.120.720:FF:000011, FunFam:1.10.10.820:FF:000001,
                FunFam:1.20.5.190:FF:000018, FunFam:1.20.5.190:FF:000001,
                Gene3D:1.10.10.820, Gene3D:1.20.5.190, Gene3D:1.20.58.530,
                Gene3D:6.20.240.20, Gene3D:3.40.850.10, Gene3D:1.20.120.720,
                InterPro:IPR002710, InterPro:IPR000048, InterPro:IPR036961,
                InterPro:IPR001609, InterPro:IPR004009, InterPro:IPR037975,
                InterPro:IPR036018, InterPro:IPR027417, PANTHER:PTHR13140,
                PANTHER:PTHR13140:SF836, Pfam:PF01843, Pfam:PF00612, Pfam:PF00063,
                Pfam:PF02736, PRINTS:PR00193, SMART:SM01132, SMART:SM00015,
                SMART:SM00242, SUPFAM:SSF52540, PROSITE:PS51126, PROSITE:PS50096,
                PROSITE:PS51456, PROSITE:PS51844
    KEYWORDS    3D-structure; Actin-binding; Alternative splicing; ATP-binding;
                Calmodulin-binding; Coiled coil; Cytoplasm; Motor protein; Myosin;
                Nucleotide-binding; Reference proteome; Repeat.
    SOURCE      Arabidopsis thaliana (thale cress)
      ORGANISM  Arabidopsis thaliana
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; rosids; malvids; Brassicales; Brassicaceae;
                Camelineae; Arabidopsis.
    REFERENCE   1  (residues 1 to 1505)
      AUTHORS   Kinkema,M., Wang,H. and Schiefelbein,J.
      TITLE     Molecular analysis of the myosin gene family in Arabidopsis
                thaliana
      JOURNAL   Plant Mol Biol 26 (4), 1139-1153 (1994)
       PUBMED   7811972
      REMARK    NUCLEOTIDE SEQUENCE [MRNA], AND TISSUE SPECIFICITY.;
                STRAIN=cv. Columbia; TISSUE=Seedling
    REFERENCE   2  (residues 1 to 1505)
      AUTHORS   Kaneko,T., Katoh,T., Asamizu,E., Sato,S., Nakamura,Y., Kotani,H.
                and Tabata,S.
      TITLE     Direct Submission
      JOURNAL   Submitted (??-JUL-1999) to the EMBL/GenBank/DDBJ databases
      REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].;
                STRAIN=cv. Columbia
    REFERENCE   3  (residues 1 to 1505)
      AUTHORS   Cheng,C.Y., Krishnakumar,V., Chan,A.P., Thibaud-Nissen,F.,
                Schobel,S. and Town,C.D.
      TITLE     Araport11: a complete reannotation of the Arabidopsis thaliana
                reference genome
      JOURNAL   Plant J 89 (4), 789-804 (2017)
       PUBMED   27862469
      REMARK    GENOME REANNOTATION.;
                STRAIN=cv. Columbia
    REFERENCE   4  (residues 1 to 1505)
      AUTHORS   Totoki,Y., Seki,M., Ishida,J., Nakajima,M., Enju,A., Kamiya,A.,
                Narusaka,M., Shin-i,T., Nakagawa,M., Sakamoto,N., Oishi,K.,
                Kohara,Y., Kobayashi,M., Toyoda,A., Sakaki,Y., Sakurai,T., Iida,K.,
                Akiyama,K., Satou,M., Toyoda,T., Konagaya,A., Carninci,P.,
                Kawai,J., Hayashizaki,Y. and Shinozaki,K.
      TITLE     Direct Submission
      JOURNAL   Submitted (??-JUL-2006) to the EMBL/GenBank/DDBJ databases
      REMARK    NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] OF 408-1505.;
                STRAIN=cv. Columbia
    REFERENCE   5  (residues 1 to 1505)
      AUTHORS   Hodge,T. and Cope,M.J.
      TITLE     A myosin family tree
      JOURNAL   J Cell Sci 113 Pt 19, 3353-3354 (2000)
       PUBMED   10984423
      REMARK    GENE FAMILY.
    REFERENCE   6  (residues 1 to 1505)
      AUTHORS   Reddy,A.S. and Day,I.S.
      TITLE     Analysis of the myosins encoded in the recently completed
                Arabidopsis thaliana genome sequence
      JOURNAL   Genome Biol 2 (7), RESEARCH0024 (2001)
       PUBMED   11516337
      REMARK    GENE FAMILY.
    REFERENCE   7  (residues 1 to 1505)
      AUTHORS   Hashimoto,K., Igarashi,H., Mano,S., Nishimura,M., Shimmen,T. and
                Yokota,E.
      TITLE     Peroxisomal localization of a myosin XI isoform in Arabidopsis
                thaliana
      JOURNAL   Plant Cell Physiol 46 (5), 782-789 (2005)
       PUBMED   15792961
      REMARK    SUBCELLULAR LOCATION, FUNCTION, AND DISRUPTION PHENOTYPE.
    REFERENCE   8  (residues 1 to 1505)
      AUTHORS   Reisen,D. and Hanson,M.R.
      TITLE     Association of six YFP-myosin XI-tail fusions with mobile plant
                cell organelles
      JOURNAL   BMC Plant Biol 7, 6 (2007)
       PUBMED   17288617
      REMARK    SUBCELLULAR LOCATION.
                Publication Status: Online-Only
    REFERENCE   9  (residues 1 to 1505)
      AUTHORS   Li,J.F. and Nebenfuhr,A.
      TITLE     Organelle targeting of myosin XI is mediated by two globular tail
                subdomains with separate cargo binding sites
      JOURNAL   J Biol Chem 282 (28), 20593-20602 (2007)
       PUBMED   17500056
      REMARK    DOMAIN, AND SUBCELLULAR LOCATION.
    REFERENCE   10 (residues 1 to 1505)
      AUTHORS   Walter,N. and Holweg,C.L.
      TITLE     Head-neck domain of Arabidopsis myosin XI, MYA2, fused with GFP
                produces F-actin patterns that coincide with fast organelle
                streaming in different plant cells
      JOURNAL   BMC Plant Biol 8, 74 (2008)
       PUBMED   18598361
      REMARK    DOMAIN, AND ACTIN-BINDING.
                Publication Status: Online-Only
    REFERENCE   11 (residues 1 to 1505)
      AUTHORS   Hashimoto,K., Igarashi,H., Mano,S., Takenaka,C., Shiina,T.,
                Yamaguchi,M., Demura,T., Nishimura,M., Shimmen,T. and Yokota,E.
      TITLE     An isoform of Arabidopsis myosin XI interacts with small GTPases in
                its C-terminal tail region
      JOURNAL   J Exp Bot 59 (13), 3523-3531 (2008)
       PUBMED   18703495
      REMARK    INTERACTION WITH RABD1 AND RABC2A.
    REFERENCE   12 (residues 1 to 1505)
      AUTHORS   Peremyslov,V.V., Prokhnevsky,A.I., Avisar,D. and Dolja,V.V.
      TITLE     Two class XI myosins function in organelle trafficking and root
                hair development in Arabidopsis
      JOURNAL   Plant Physiol 146 (3), 1109-1116 (2008)
       PUBMED   18178669
      REMARK    DISRUPTION PHENOTYPE, AND FUNCTION.
    REFERENCE   13 (residues 1 to 1505)
      AUTHORS   Prokhnevsky,A.I., Peremyslov,V.V. and Dolja,V.V.
      TITLE     Overlapping functions of the four class XI myosins in Arabidopsis
                growth, root hair elongation, and organelle motility
      JOURNAL   Proc Natl Acad Sci U S A 105 (50), 19744-19749 (2008)
       PUBMED   19060218
      REMARK    DISRUPTION PHENOTYPE, AND FUNCTION.
    REFERENCE   14 (residues 1 to 1505)
      AUTHORS   Avisar,D., Abu-Abied,M., Belausov,E., Sadot,E., Hawes,C. and
                Sparkes,I.A.
      TITLE     A comparative study of the involvement of 17 Arabidopsis myosin
                family members on the motility of Golgi and other organelles
      JOURNAL   Plant Physiol 150 (2), 700-709 (2009)
       PUBMED   19369591
      REMARK    FUNCTION.
    REFERENCE   15 (residues 1 to 1505)
      AUTHORS   Peremyslov,V.V., Prokhnevsky,A.I. and Dolja,V.V.
      TITLE     Class XI myosins are required for development, cell expansion, and
                F-Actin organization in Arabidopsis
      JOURNAL   Plant Cell 22 (6), 1883-1897 (2010)
       PUBMED   20581304
      REMARK    FUNCTION.
    REFERENCE   16 (residues 1 to 1505)
      AUTHORS   Peremyslov,V.V., Mockler,T.C., Filichkin,S.A., Fox,S.E.,
                Jaiswal,P., Makarova,K.S., Koonin,E.V. and Dolja,V.V.
      TITLE     Expression, splicing, and evolution of the myosin gene family in
                plants
      JOURNAL   Plant Physiol 155 (3), 1191-1204 (2011)
       PUBMED   21233331
      REMARK    GENE FAMILY, AND NOMENCLATURE.
    REFERENCE   17 (residues 1 to 1505)
      AUTHORS   Ojangu,E.L., Tanner,K., Pata,P., Jarve,K., Holweg,C.L., Truve,E.
                and Paves,H.
      TITLE     Myosins XI-K, XI-1, and XI-2 are required for development of
                pavement cells, trichomes, and stigmatic papillae in Arabidopsis
      JOURNAL   BMC Plant Biol 12, 81 (2012)
       PUBMED   22672737
      REMARK    FUNCTION.
                Publication Status: Online-Only
    REFERENCE   18 (residues 1 to 1505)
      AUTHORS   Avisar,D., Abu-Abied,M., Belausov,E. and Sadot,E.
      TITLE     Myosin XIK is a major player in cytoplasm dynamics and is regulated
                by two amino acids in its tail
      JOURNAL   J Exp Bot 63 (1), 241-249 (2012)
       PUBMED   21914656
      REMARK    FUNCTION.
    COMMENT     On or before Jul 2, 2013 this sequence version replaced
                gi:122229993, gi:75101940.
                [FUNCTION] Myosin heavy chain that is required for the cell
                cycle-regulated transport of various organelles and proteins for
                their segregation. Functions by binding with its tail domain to
                receptor proteins on organelles and exerting force with its
                N-terminal motor domain against actin filaments, thereby
                transporting its cargo along polarized actin cables. Involved in
                the tip growth of root hair cells. Plays a major role in
                trafficking of Golgi stacks, mitochondria and peroxisomes during
                root hair development. Targets the peroxisome through an
                interaction with RABC2A. Required for development of pavement
                cells, trichomes, and stigmatic papillae.
                {ECO:0000269|PubMed:15792961, ECO:0000269|PubMed:18178669,
                ECO:0000269|PubMed:19060218, ECO:0000269|PubMed:19369591,
                ECO:0000269|PubMed:20581304, ECO:0000269|PubMed:21914656,
                ECO:0000269|PubMed:22672737}.
                [SUBUNIT] Homodimer (By similarity). Interacts with RABC2A and
                RABD1. {ECO:0000250, ECO:0000269|PubMed:18703495}.
                [INTERACTION] Q9LKB9; O49841: RABC2A; NbExp=3; IntAct=EBI-2009528,
                EBI-2009559; Q9LKB9; Q9ZRE2: RABD1; NbExp=3; IntAct=EBI-2009528,
                EBI-2009542.
                [SUBCELLULAR LOCATION] Cytoplasm {ECO:0000269|PubMed:15792961,
                ECO:0000269|PubMed:17288617, ECO:0000269|PubMed:17500056}.
                Note=Colocalizes with peroxisome and cytoplasmic vesicles.
                [ALTERNATIVE PRODUCTS] Event=Alternative splicing; Named
                isoforms=1; Comment=A number of isoforms are produced. According to
                EST sequences.; Name=1; IsoId=Q9LKB9-1; Sequence=Displayed.
                [TISSUE SPECIFICITY] Expressed in flowers, leaves, roots and stems.
                {ECO:0000269|PubMed:7811972}.
                [DOMAIN] Head-neck domain associates with cytoplasmic
                (transvacuolar) F-actin in areas coinciding with the tracks of fast
                organelles. {ECO:0000269|PubMed:17500056,
                ECO:0000269|PubMed:18598361}.
                [DOMAIN] IQ domain mediates interaction with calmodulin.
                {ECO:0000250}.
                [DOMAIN] The tail domain is a globular cargo-binding domain.
                {ECO:0000250}.
                [DISRUPTION PHENOTYPE] Impaired growth of root hair cells
                (PubMed:18178669). No visible phenotype (PubMed:15792961).
                {ECO:0000269|PubMed:15792961, ECO:0000269|PubMed:18178669,
                ECO:0000269|PubMed:19060218}.
                [SIMILARITY] Belongs to the TRAFAC class myosin-kinesin ATPase
                superfamily. Myosin family. Plant myosin class XI subfamily.
                {ECO:0000305}.
                [SEQUENCE CAUTION] Sequence=CAA84066.1; Type=Frameshift;
                Evidence={ECO:0000305}.
    FEATURES             Location/Qualifiers
         source          1..1505
                         /organism="Arabidopsis thaliana"
                         /db_xref="taxon:3702"
         gene            1..1505
                         /gene="XI-2"
                         /locus_tag="At5g43900"
                         /gene_synonym="F6B6.4"
                         /gene_synonym="MYA2"
         Protein         1..1505
                         /product="Myosin-6"
                         /note="AtMYA2"
                         /UniProtKB_evidence="Evidence at protein level"
         Region          1..1505
                         /region_name="Mature chain"
                         /note="Myosin-6. /id=PRO_0000422861."
         Region          8..57
                         /region_name="Domain"
                         /note="Myosin N-terminal SH3-like.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU01190."
         Region          8..48
                         /region_name="Myosin_N"
                         /note="Myosin N-terminal SH3-like domain; pfam02736"
                         /db_xref="CDD:460670"
         Region          12..16
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          18..30
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          32..39
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          44..48
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          49..51
                         /region_name="Hydrogen bonded turn"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          62..731
                         /region_name="Domain"
                         /note="Myosin motor.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00782."
         Region          67..69
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          75..87
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          76..719
                         /region_name="MYSc_Myo11"
                         /note="class XI myosin, motor domain; cd01384"
                         /db_xref="CDD:276835"
         Region          92..95
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          98..102
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            order(103..112,156..163,204..214,434..439)
                         /site_type="other"
                         /note="ATP binding site [chemical binding]"
                         /db_xref="CDD:276835"
         Site            103..112
                         /site_type="other"
                         /note="purine-binding loop"
                         /db_xref="CDD:276835"
         Region          110..112
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          114..119
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          125..127
                         /region_name="Hydrogen bonded turn"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          132..146
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          136
                         /region_name="Conflict"
                         /note="A -> P (in Ref. 1; CAA84066).
                         /evidence=ECO:0000305."
         Region          150..155
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            156..163
                         /site_type="other"
                         /note="P-loop"
                         /db_xref="CDD:276835"
         Region          162..177
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          187..202
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            204..214
                         /site_type="other"
                         /note="switch I region"
                         /db_xref="CDD:276835"
         Region          212..223
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          229..237
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          242..244
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          249..251
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          255..261
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          265..267
                         /region_name="Hydrogen bonded turn"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          269..272
                         /region_name="Hydrogen bonded turn"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          276..278
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          280..283
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          295..309
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          313..330
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          334..336
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          343..345
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          345
                         /region_name="Conflict"
                         /note="P -> S (in Ref. 1; CAA84066).
                         /evidence=ECO:0000305."
         Region          348..360
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          365..373
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          390..420
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          427..434
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            434..439
                         /site_type="other"
                         /note="switch II region"
                         /db_xref="CDD:276835"
         Region          446..466
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            order(463..474,477..486)
                         /site_type="active"
                         /note="relay loop [active]"
                         /db_xref="CDD:276835"
         Region          468..476
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          490..497
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          495..529
                         /region_name="Region of interest in the sequence"
                         /note="Actin-binding. /evidence=ECO:0000255."
         Region          499..501
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          503..512
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          513..515
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          518..528
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          529..531
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          531..554
                         /region_name="Region of interest in the sequence"
                         /note="Actin-binding. /evidence=ECO:0000255."
         Region          533..536
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          542..549
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          552..557
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          558..560
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          561..565
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          571..578
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          583..586
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          589..612
                         /region_name="Region of interest in the sequence"
                         /note="Actin-binding. /evidence=ECO:0000255."
         Region          593..595
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          604..619
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          612..634
                         /region_name="Region of interest in the sequence"
                         /note="Actin-binding. /evidence=ECO:0000250."
         Region          622..630
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          643..652
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            654..663
                         /site_type="other"
                         /note="SH1 helix"
                         /db_xref="CDD:276835"
         Region          655..664
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Site            order(666..685,708..719)
                         /site_type="other"
                         /note="converter subdomain"
                         /db_xref="CDD:276835"
         Region          668..671
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          672..679
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          680..682
                         /region_name="Hydrogen bonded turn"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          694..704
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          713..718
                         /region_name="Beta-strand region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          722..735
                         /region_name="Helical region"
                         /note="/evidence=ECO:0007829|PDB:7DHW."
         Region          734..763
                         /region_name="Domain"
                         /note="IQ 1.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00116."
         Region          757..786
                         /region_name="Domain"
                         /note="IQ 2.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00116."
         Region          782..811
                         /region_name="Domain"
                         /note="IQ 3.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00116."
         Region          805..834
                         /region_name="Domain"
                         /note="IQ 4.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00116."
         Region          829..851
                         /region_name="IQ"
                         /note="Calmodulin-binding motif; smart00015"
                         /db_xref="CDD:197470"
         Region          830..859
                         /region_name="Domain"
                         /note="IQ 5.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00116."
         Region          853..882
                         /region_name="Domain"
                         /note="IQ 6.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00116."
         Region          <875..>1051
                         /region_name="Myosin_tail_1"
                         /note="Myosin tail; pfam01576"
                         /db_xref="CDD:460256"
         Region          883..1048
                         /region_name="Coiled-coil region"
                         /note="/evidence=ECO:0000255."
         Region          1099..1475
                         /region_name="MyosinXI_CBD"
                         /note="cargo binding domain of myosin XI; cd15475"
                         /db_xref="CDD:271259"
         Region          1148..1452
                         /region_name="Domain"
                         /note="Dilute.
                         /evidence=ECO:0000255|PROSITE-ProRule:PRU00503."
         Region          1442
                         /region_name="Conflict"
                         /note="V -> A (in Ref. 4; BAF00999).
                         /evidence=ECO:0000305."
    ORIGIN      
            1 mvanfnpsvg sfvwvedpde awidgevvqv ngdeikvlct sgkhvvtkis naypkdveap
           61 asgvddmtrl aylhepgvlq nlhsrydine iytytgsili avnpfrrlph lysshmmaqy
          121 kgaslgelsp hpfavadaay rqmindgvsq silvsgesga gktestkllm rylaymggra
          181 aaegrsveqk vlesnpvlea fgnaktvrnn nssrfgkfve iqfdekgris gaairtylle
          241 rsrvcqvsdp ernyhcfyml caapqedvkk fkleepkkyh ylnqskclel dsindaeeyh
          301 atrramdvvg isteeqdaif svvaailhig niefakgeei dssipkddks lfhlktaael
          361 lscdekaled slckrimvtr detitktldp eaatlsrdal akvmysrlfd wlvdkinssi
          421 gqdhdskyli gvldiygfes fktnsfeqfc inltneklqq hfnqhvfkme qeeykkeein
          481 wsyiefvdnq dildliekkp ggiialldea cmfprsthet faqklyqtfk thkrftkpkl
          541 arsdftichy agdvtyqtel fldknkdyvi aehqallnss scsfvaslfp pmsddskqsk
          601 fssigtrfkq qlvslleiln ttephyirci kpnnllkpgi fenenilqql rcggvmeair
          661 iscagyptrk hfdeflarfg ilapevlvkn sddpaackkl ldkvglegyq igktkvflra
          721 gqmadldtrr tevlgrsasi iqrkvrsyla kksfivlrns akqiqsvcrg ylarsvyegm
          781 rreaaalkiq rdlrrflark aytelysaav svqagmrgmv arkelcfrrq tkaaiiiqtw
          841 crgylarlhy rklkkaaitt qcawrskvar gelrklkmaa retgalqaak nklekqveel
          901 twrlqlekri rtdleeakkq esakaqssle elqlkckete allikereaa kkiaetapii
          961 keipvvdqel mdkitnenek lksmvsslem kigetekklq ettkisqdrl nqaleaeskl
         1021 vklktamqrl eekildmeae kkimhqqtis tpvrtnlghp ptapvknlen ghqtnlekef
         1081 neaefttpvd gkagksaaer qimnvdalid cvkdnigfsn gkpvaaftiy kcllhwkcfe
         1141 sektnvfdrl iqmigsaien eddnshlayw ltstsallfl lqkslktngs gatqskkppa
         1201 stslfgrmam sfrsspasgn laaaaeaaal avvrpveaky pallfkqqla ayvekmfgmv
         1261 rdnlkrelst llslciqapr sskggmlrsg rsfgkdspav hwqsiidgln sllvtlkenh
         1321 vplvliqkiy sqtfsyinvq lfnslllrke cctfsngefv ksglaelelw ccqakeysgp
         1381 sweelkhirq avgflvihqk yrisydeian dlcpvlsvqq lyrictlywd dsyntrsvsq
         1441 evissmrtlm teesndadsd sflldddssi pfsiddisss meekdfvgik paeellenpa
         1501 fvflh
    //
    
    
    Fetching domains for: XP_016503116.1
    LOCUS       XP_016503116             202 aa            linear   PLN 09-DEC-2024
    DEFINITION  ras-related protein RABD1 [Nicotiana tabacum].
    ACCESSION   XP_016503116
    VERSION     XP_016503116.1
    DBLINK      BioProject: PRJNA319578
    DBSOURCE    REFSEQ: accession XM_016647630.1
    KEYWORDS    RefSeq.
    SOURCE      Nicotiana tabacum (common tobacco)
      ORGANISM  Nicotiana tabacum
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Solanales; Solanaceae;
                Nicotianoideae; Nicotianeae; Nicotiana.
    COMMENT     MODEL REFSEQ:  This record is predicted by automated computational
                analysis. This record is derived from a genomic sequence
                (NW_015827563.1) annotated using gene prediction method: Gnomon,
                supported by EST evidence.
                Also see:
                    Documentation of NCBI's Annotation Process
                
                ##Genome-Annotation-Data-START##
                Annotation Provider         :: NCBI RefSeq
                Annotation Status           :: Updated annotation
                Annotation Name             :: GCF_000715135.1-RS_2024_11
                Annotation Pipeline         :: NCBI eukaryotic genome annotation
                                               pipeline
                Annotation Software Version :: 10.3
                Annotation Method           :: Best-placed RefSeq; Gnomon;
                                               tRNAscan-SE
                Features Annotated          :: Gene; mRNA; CDS; ncRNA
                Annotation Date             :: 11/26/2024
                ##Genome-Annotation-Data-END##
                COMPLETENESS: full length.
    FEATURES             Location/Qualifiers
         source          1..202
                         /organism="Nicotiana tabacum"
                         /cultivar="TN90"
                         /db_xref="taxon:4097"
                         /chromosome="Unknown"
         Protein         1..202
                         /product="ras-related protein RABD1"
                         /calculated_mol_wt=22316
         Region          7..172
                         /region_name="Rab1_Ypt1"
                         /note="Rab GTPase family 1 includes the yeast homolog
                         Ypt1; cd01869"
                         /db_xref="CDD:206661"
         Site            7..10
                         /site_type="other"
                         /note="Rab subfamily motif 1 (RabSF1)"
                         /db_xref="CDD:206661"
         Site            15..22
                         /site_type="other"
                         /note="G1 box"
                         /db_xref="CDD:206661"
         Site            order(16,18..23,33,121..122,124,151..153)
                         /site_type="other"
                         /note="GTP/Mg2+ binding site [chemical binding]"
                         /db_xref="CDD:206661"
         Site            order(23..34,37..39)
                         /site_type="other"
                         /note="Rab subfamily motif 2 (RabSF2)"
                         /db_xref="CDD:206661"
         Site            order(33,38..46)
                         /site_type="other"
                         /note="Switch I region"
                         /db_xref="CDD:206661"
         Site            order(38,42..49,56,58)
                         /site_type="other"
                         /note="putative GEF interaction site [polypeptide
                         binding]"
                         /db_xref="CDD:206661"
         Site            40
                         /site_type="other"
                         /note="G2 box"
                         /db_xref="CDD:206661"
         Site            order(41,43..45,60,62,69..70,73,77,79..82,169)
                         /site_type="active"
                         /note="putative effector interaction site [active]"
                         /db_xref="CDD:206661"
         Site            order(41..42,44,62..63,70,72,74..76,79)
                         /site_type="active"
                         /note="GDI interaction site [active]"
                         /db_xref="CDD:206661"
         Site            41..45
                         /site_type="other"
                         /note="Rab family motif 1 (RabF1)"
                         /db_xref="CDD:206661"
         Site            58..62
                         /site_type="other"
                         /note="Rab family motif 2 (RabF2)"
                         /db_xref="CDD:206661"
         Site            63..66
                         /site_type="other"
                         /note="G3 box"
                         /db_xref="CDD:206661"
         Site            order(66,68..78)
                         /site_type="other"
                         /note="Switch II region"
                         /db_xref="CDD:206661"
         Site            69..74
                         /site_type="other"
                         /note="Rab family motif 3 (RabF3)"
                         /db_xref="CDD:206661"
         Site            77..81
                         /site_type="other"
                         /note="Rab family motif 4 (RabF4)"
                         /db_xref="CDD:206661"
         Site            86..91
                         /site_type="other"
                         /note="Rab family motif 5 (RabF5)"
                         /db_xref="CDD:206661"
         Site            112..118
                         /site_type="other"
                         /note="Rab subfamily motif 3 (RabSF3)"
                         /db_xref="CDD:206661"
         Site            121..124
                         /site_type="other"
                         /note="G4 box"
                         /db_xref="CDD:206661"
         Site            151..153
                         /site_type="other"
                         /note="G5 box"
                         /db_xref="CDD:206661"
         Site            167..172
                         /site_type="other"
                         /note="Rab subfamily motif 4 (RabSF4)"
                         /db_xref="CDD:206661"
         CDS             1..202
                         /gene="LOC107821208"
                         /coded_by="XM_016647630.1:127..735"
                         /db_xref="GeneID:107821208"
    ORIGIN      
            1 msneydylfk llligdssvg ksclllrfad dsyvdsyist igvdfkirtv eldgktiklq
           61 iwdtagqerf rtitssyyrg ahgiiivydv temesfnnvk qwlneidrya nesvckllvg
          121 nkcdlvenkv vdtqtgkala delgipflet sakdsinveq afltmageik kkmgnqpaga
          181 kksgstvqik gqpieqksnc cg
    //
    
    
    Fetching domains for: XP_016446075.1
    LOCUS       XP_016446075             202 aa            linear   PLN 09-DEC-2024
    DEFINITION  ras-related protein RABD1 [Nicotiana tabacum].
    ACCESSION   XP_016446075
    VERSION     XP_016446075.1
    DBLINK      BioProject: PRJNA319578
    DBSOURCE    REFSEQ: accession XM_016590589.1
    KEYWORDS    RefSeq.
    SOURCE      Nicotiana tabacum (common tobacco)
      ORGANISM  Nicotiana tabacum
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Solanales; Solanaceae;
                Nicotianoideae; Nicotianeae; Nicotiana.
    COMMENT     MODEL REFSEQ:  This record is predicted by automated computational
                analysis. This record is derived from a genomic sequence
                (NW_015891302.1) annotated using gene prediction method: Gnomon,
                supported by EST evidence.
                Also see:
                    Documentation of NCBI's Annotation Process
                
                ##Genome-Annotation-Data-START##
                Annotation Provider         :: NCBI RefSeq
                Annotation Status           :: Updated annotation
                Annotation Name             :: GCF_000715135.1-RS_2024_11
                Annotation Pipeline         :: NCBI eukaryotic genome annotation
                                               pipeline
                Annotation Software Version :: 10.3
                Annotation Method           :: Best-placed RefSeq; Gnomon;
                                               tRNAscan-SE
                Features Annotated          :: Gene; mRNA; CDS; ncRNA
                Annotation Date             :: 11/26/2024
                ##Genome-Annotation-Data-END##
                COMPLETENESS: full length.
    FEATURES             Location/Qualifiers
         source          1..202
                         /organism="Nicotiana tabacum"
                         /cultivar="TN90"
                         /db_xref="taxon:4097"
                         /chromosome="Unknown"
         Protein         1..202
                         /product="ras-related protein RABD1"
                         /calculated_mol_wt=22330
         Region          7..172
                         /region_name="Rab1_Ypt1"
                         /note="Rab GTPase family 1 includes the yeast homolog
                         Ypt1; cd01869"
                         /db_xref="CDD:206661"
         Site            7..10
                         /site_type="other"
                         /note="Rab subfamily motif 1 (RabSF1)"
                         /db_xref="CDD:206661"
         Site            15..22
                         /site_type="other"
                         /note="G1 box"
                         /db_xref="CDD:206661"
         Site            order(16,18..23,33,121..122,124,151..153)
                         /site_type="other"
                         /note="GTP/Mg2+ binding site [chemical binding]"
                         /db_xref="CDD:206661"
         Site            order(23..34,37..39)
                         /site_type="other"
                         /note="Rab subfamily motif 2 (RabSF2)"
                         /db_xref="CDD:206661"
         Site            order(33,38..46)
                         /site_type="other"
                         /note="Switch I region"
                         /db_xref="CDD:206661"
         Site            order(38,42..49,56,58)
                         /site_type="other"
                         /note="putative GEF interaction site [polypeptide
                         binding]"
                         /db_xref="CDD:206661"
         Site            40
                         /site_type="other"
                         /note="G2 box"
                         /db_xref="CDD:206661"
         Site            order(41,43..45,60,62,69..70,73,77,79..82,169)
                         /site_type="active"
                         /note="putative effector interaction site [active]"
                         /db_xref="CDD:206661"
         Site            order(41..42,44,62..63,70,72,74..76,79)
                         /site_type="active"
                         /note="GDI interaction site [active]"
                         /db_xref="CDD:206661"
         Site            41..45
                         /site_type="other"
                         /note="Rab family motif 1 (RabF1)"
                         /db_xref="CDD:206661"
         Site            58..62
                         /site_type="other"
                         /note="Rab family motif 2 (RabF2)"
                         /db_xref="CDD:206661"
         Site            63..66
                         /site_type="other"
                         /note="G3 box"
                         /db_xref="CDD:206661"
         Site            order(66,68..78)
                         /site_type="other"
                         /note="Switch II region"
                         /db_xref="CDD:206661"
         Site            69..74
                         /site_type="other"
                         /note="Rab family motif 3 (RabF3)"
                         /db_xref="CDD:206661"
         Site            77..81
                         /site_type="other"
                         /note="Rab family motif 4 (RabF4)"
                         /db_xref="CDD:206661"
         Site            86..91
                         /site_type="other"
                         /note="Rab family motif 5 (RabF5)"
                         /db_xref="CDD:206661"
         Site            112..118
                         /site_type="other"
                         /note="Rab subfamily motif 3 (RabSF3)"
                         /db_xref="CDD:206661"
         Site            121..124
                         /site_type="other"
                         /note="G4 box"
                         /db_xref="CDD:206661"
         Site            151..153
                         /site_type="other"
                         /note="G5 box"
                         /db_xref="CDD:206661"
         Site            167..172
                         /site_type="other"
                         /note="Rab subfamily motif 4 (RabSF4)"
                         /db_xref="CDD:206661"
         CDS             1..202
                         /gene="LOC107771246"
                         /coded_by="XM_016590589.1:97..705"
                         /db_xref="GeneID:107771246"
    ORIGIN      
            1 msneydylfk llligdssvg ksclllrfad dsyvdsyist igvdfkirtv eldgktiklq
           61 iwdtagqerf rtitssyyrg ahgiiivydv temesfnnvk qwlneidrya nesvckllvg
          121 nkcdlvenkv vdtqtgkala delgipflet sakdsinveq afltmageik kkmgnqpaga
          181 kktgstvqik gqpieqksnc cg
    //
    
    
    Fetching domains for: XP_016461491.1
    LOCUS       XP_016461491             202 aa            linear   PLN 09-DEC-2024
    DEFINITION  ras-related protein RABD1 [Nicotiana tabacum].
    ACCESSION   XP_016461491
    VERSION     XP_016461491.1
    DBLINK      BioProject: PRJNA319578
    DBSOURCE    REFSEQ: accession XM_016606005.1
    KEYWORDS    RefSeq.
    SOURCE      Nicotiana tabacum (common tobacco)
      ORGANISM  Nicotiana tabacum
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Solanales; Solanaceae;
                Nicotianoideae; Nicotianeae; Nicotiana.
    COMMENT     MODEL REFSEQ:  This record is predicted by automated computational
                analysis. This record is derived from a genomic sequence
                (NW_015788348.1) annotated using gene prediction method: Gnomon,
                supported by EST evidence.
                Also see:
                    Documentation of NCBI's Annotation Process
                
                ##Genome-Annotation-Data-START##
                Annotation Provider         :: NCBI RefSeq
                Annotation Status           :: Updated annotation
                Annotation Name             :: GCF_000715135.1-RS_2024_11
                Annotation Pipeline         :: NCBI eukaryotic genome annotation
                                               pipeline
                Annotation Software Version :: 10.3
                Annotation Method           :: Best-placed RefSeq; Gnomon;
                                               tRNAscan-SE
                Features Annotated          :: Gene; mRNA; CDS; ncRNA
                Annotation Date             :: 11/26/2024
                ##Genome-Annotation-Data-END##
                COMPLETENESS: full length.
    FEATURES             Location/Qualifiers
         source          1..202
                         /organism="Nicotiana tabacum"
                         /cultivar="TN90"
                         /db_xref="taxon:4097"
                         /chromosome="Unknown"
         Protein         1..202
                         /product="ras-related protein RABD1"
                         /calculated_mol_wt=22348
         Region          7..172
                         /region_name="Rab1_Ypt1"
                         /note="Rab GTPase family 1 includes the yeast homolog
                         Ypt1; cd01869"
                         /db_xref="CDD:206661"
         Site            7..10
                         /site_type="other"
                         /note="Rab subfamily motif 1 (RabSF1)"
                         /db_xref="CDD:206661"
         Site            15..22
                         /site_type="other"
                         /note="G1 box"
                         /db_xref="CDD:206661"
         Site            order(16,18..23,33,121..122,124,151..153)
                         /site_type="other"
                         /note="GTP/Mg2+ binding site [chemical binding]"
                         /db_xref="CDD:206661"
         Site            order(23..34,37..39)
                         /site_type="other"
                         /note="Rab subfamily motif 2 (RabSF2)"
                         /db_xref="CDD:206661"
         Site            order(33,38..46)
                         /site_type="other"
                         /note="Switch I region"
                         /db_xref="CDD:206661"
         Site            order(38,42..49,56,58)
                         /site_type="other"
                         /note="putative GEF interaction site [polypeptide
                         binding]"
                         /db_xref="CDD:206661"
         Site            40
                         /site_type="other"
                         /note="G2 box"
                         /db_xref="CDD:206661"
         Site            order(41,43..45,60,62,69..70,73,77,79..82,169)
                         /site_type="active"
                         /note="putative effector interaction site [active]"
                         /db_xref="CDD:206661"
         Site            order(41..42,44,62..63,70,72,74..76,79)
                         /site_type="active"
                         /note="GDI interaction site [active]"
                         /db_xref="CDD:206661"
         Site            41..45
                         /site_type="other"
                         /note="Rab family motif 1 (RabF1)"
                         /db_xref="CDD:206661"
         Site            58..62
                         /site_type="other"
                         /note="Rab family motif 2 (RabF2)"
                         /db_xref="CDD:206661"
         Site            63..66
                         /site_type="other"
                         /note="G3 box"
                         /db_xref="CDD:206661"
         Site            order(66,68..78)
                         /site_type="other"
                         /note="Switch II region"
                         /db_xref="CDD:206661"
         Site            69..74
                         /site_type="other"
                         /note="Rab family motif 3 (RabF3)"
                         /db_xref="CDD:206661"
         Site            77..81
                         /site_type="other"
                         /note="Rab family motif 4 (RabF4)"
                         /db_xref="CDD:206661"
         Site            86..91
                         /site_type="other"
                         /note="Rab family motif 5 (RabF5)"
                         /db_xref="CDD:206661"
         Site            112..118
                         /site_type="other"
                         /note="Rab subfamily motif 3 (RabSF3)"
                         /db_xref="CDD:206661"
         Site            121..124
                         /site_type="other"
                         /note="G4 box"
                         /db_xref="CDD:206661"
         Site            151..153
                         /site_type="other"
                         /note="G5 box"
                         /db_xref="CDD:206661"
         Site            167..172
                         /site_type="other"
                         /note="Rab subfamily motif 4 (RabSF4)"
                         /db_xref="CDD:206661"
         CDS             1..202
                         /gene="LOC107784791"
                         /coded_by="XM_016606005.1:146..754"
                         /db_xref="GeneID:107784791"
    ORIGIN      
            1 msneydylfk llligdssvg ksclllrfad dsyvdsyist igvdfkirtv eldgktiklq
           61 iwdtagqerf rtitssyyrg ahgiiivydv tekesfdnvk qwlseidrya nesvckllvg
          121 nkcdlvenkv vdtqtakafa delgipflet sakdsinveq afltmageik kkmgsqpagt
          181 knsgrgvqik gqpieqksnc cg
    //
    
    
    Fetching domains for: XP_070037692.1
    LOCUS       XP_070037692             202 aa            linear   PLN 05-DEC-2024
    DEFINITION  ras-related protein RABD1-like [Nicotiana tomentosiformis].
    ACCESSION   XP_070037692
    VERSION     XP_070037692.1
    DBLINK      BioProject: PRJNA257218
    DBSOURCE    REFSEQ: accession XM_070181591.1
    KEYWORDS    RefSeq.
    SOURCE      Nicotiana tomentosiformis
      ORGANISM  Nicotiana tomentosiformis
                Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
                Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae;
                Pentapetalae; asterids; lamiids; Solanales; Solanaceae;
                Nicotianoideae; Nicotianeae; Nicotiana.
    COMMENT     MODEL REFSEQ:  This record is predicted by automated computational
                analysis. This record is derived from a genomic sequence
                (NC_090818) annotated using gene prediction method: Gnomon.
                Also see:
                    Documentation of NCBI's Annotation Process
                
                ##Genome-Annotation-Data-START##
                Annotation Provider         :: NCBI RefSeq
                Annotation Status           :: Full annotation
                Annotation Name             :: GCF_000390325.3-RS_2024_10
                Annotation Pipeline         :: NCBI eukaryotic genome annotation
                                               pipeline
                Annotation Software Version :: 10.3
                Annotation Method           :: Best-placed RefSeq; Gnomon;
                                               cmsearch; tRNAscan-SE
                Features Annotated          :: Gene; mRNA; CDS; ncRNA
                Annotation Date             :: 10/03/2024
                ##Genome-Annotation-Data-END##
                COMPLETENESS: full length.
    FEATURES             Location/Qualifiers
         source          1..202
                         /organism="Nicotiana tomentosiformis"
                         /bio_material="USDA:TW 142"
                         /db_xref="taxon:4098"
                         /chromosome="7"
                         /tissue_type="leaf"
         Protein         1..202
                         /product="ras-related protein RABD1-like"
                         /calculated_mol_wt=22333
         CDS             1..202
                         /gene="LOC104120579"
                         /coded_by="XM_070181591.1:35..643"
                         /db_xref="GeneID:104120579"
    ORIGIN      
            1 msneydylfk llligdssvg ksclllrfad dsyvdsyist igvdfkirtv eldgktiklq
           61 iwdtagqerf rtitssyyrg ahgiiivydv tekesfdnvk qwlseidrya nesvckllvg
          121 nkcdlvenkv vdtqiakafa delgipflet sakdsinveq afltmageik kkmgnqpagt
          181 kksgsgvqik gqpieqksnc cg
    //
    
    
    Visualizing motifs...


    /opt/anaconda3/lib/python3.12/site-packages/Bio/Align/AlignInfo.py:364: BiopythonDeprecationWarning: The `pos_specific_score_matrix` method is deprecated and will be removed in a future release of Biopython. As an alternative, you can convert the multiple sequence alignment object to a new-style Alignment object by via its `.alignment` property, and then create a Motif object. For example, for a multiple sequence alignment `msa` of DNA nucleotides, you would do: 
    >>> alignment = msa.alignment
    >>> from Bio.motifs import Motif
    >>> motif = Motif('ACGT', alignment)
    >>> counts = motif.counts
    
    The `counts` object contains the same information as the PSSM returned by `pos_specific_score_matrix`, but note that the indices are reversed:
    
    >>> counts[letter][i] == pssm[index][letter]
    True
    
    If your multiple sequence alignment object was obtained using Bio.AlignIO, then you can obtain a new-style Alignment object directly by using Bio.Align.read instead of Bio.AlignIO.read, or Bio.Align.parse instead of Bio.AlignIO.parse.
      warnings.warn(
    /opt/anaconda3/lib/python3.12/site-packages/Bio/Align/AlignInfo.py:62: BiopythonDeprecationWarning: The `dumb_consensus` method is deprecated and will be removed in a future release of Biopython. As an alternative, you can convert the multiple sequence alignment object to a new-style Alignment object by via its `.alignment` property, and then create a Motif object. You can then use the `.consensus` or `.degenerate_consensus` property of the Motif object to get a consensus sequence. For more control over how the consensus sequence is calculated, you can call the `calculate_consensus` method on the `.counts` property of the Motif object. This is an example for a multiple sequence alignment `msa` of DNA nucleotides:
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Align import MultipleSeqAlignment
    >>> from Bio.Align.AlignInfo import SummaryInfo
    >>> msa = MultipleSeqAlignment([SeqRecord(Seq('ACGT')),
    ...                             SeqRecord(Seq('ATGT')),
    ...                             SeqRecord(Seq('ATGT'))])
    >>> summary = SummaryInfo(msa)
    >>> dumb_consensus = summary.dumb_consensus(ambiguous='N')
    >>> print(dumb_consensus)
    ANGT
    >>> alignment = msa.alignment
    >>> from Bio.motifs import Motif
    >>> motif = Motif('ACGT', alignment)
    >>> print(motif.consensus)
    ATGT
    >>> print(motif.degenerate_consensus)
    AYGT
    >>> counts = motif.counts
    >>> consensus = counts.calculate_consensus(identity=0.7)
    >>> print(consensus)
    ANGT
    
    If your multiple sequence alignment object was obtained using Bio.AlignIO, then you can obtain a new-style Alignment object directly by using Bio.Align.read instead of Bio.AlignIO.read, or Bio.Align.parse instead of Bio.AlignIO.parse.
      warnings.warn(
    /opt/anaconda3/lib/python3.12/site-packages/Bio/Align/AlignInfo.py:718: BiopythonDeprecationWarning: The `PSSM` class is deprecated and will be removed in a future release of Biopython. As an alternative, you can convert the multiple sequence alignment object to a new-style Alignment object by via its `.alignment` property, and then create a Motif object. For example, for a multiple sequence alignment `msa` of DNA nucleotides, you would do: 
    >>> alignment = msa.alignment
    >>> from Bio.motifs import Motif
    >>> motif = Motif('ACGT', alignment)
    >>> counts = motif.counts
    
    The `counts` object contains the same information as the PSSM returned by `pos_specific_score_matrix`, but note that the indices are reversed:
    
    >>> counts[letter][i] == pssm[index][letter]
    True
    
    If your multiple sequence alignment object was obtained using Bio.AlignIO, then you can obtain a new-style Alignment object directly by using Bio.Align.read instead of Bio.AlignIO.read, or Bio.Align.parse instead of Bio.AlignIO.parse.
      warnings.warn(



    
![png](output_22_4.png)
    


    Visualizing sequence motifs using sequence logo...



    <Figure size 1400x800 with 0 Axes>



    
![png](output_22_7.png)
    


## 3.Executing the Main Script 

**Aim:** Execute the pipeline where all functions have been saved to *main()*

**Overview of Final Pipeline**
1. **Query Sequence Input:**
   - Protein coding query sequences (**RABC2A** and **RABD1**) defined and saved in FASTA format (these can be replaced with any query sequences)
2. **Retreieve Homologous Sequences**
   - Using the ***align_sequences()*** function
3. **Sequence Alignment**
   - Homologous sequences were aligned using the ***align_sequences()*** function
4. **Constructing a Phylogenetic Tree**
   - Visualize the genetic similarity between queries and homologues using the ***build_phylogenetic_tree()*** function
5. **Conserved domain analysis**
   - Using the ***analyze_conserved_domains()*** function
6. **Visualize Motifs**
   - Motifs visualized through a PSSM heatmap using the ***visualize_motifs()*** function
7. **Sequence Logo**
   - Sequence Logo created to display the conserved features using the ***visualize_sequence_logo()*** function 


```python

```
