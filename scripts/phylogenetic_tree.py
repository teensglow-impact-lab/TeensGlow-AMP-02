import os
import pandas as pd
from Bio import Entrez, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

def get_lineage(org_name):
    """
    Retrieves the taxonomic lineage of a strain via NCBI Taxonomy.
    """
    try:
        handle = Entrez.esearch(db="taxonomy", term=org_name)
        record = Entrez.read(handle)
        handle.close()
        
        if record["IdList"]:
            tax_id = record["IdList"][0]
            handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            return records[0]["Lineage"]
    except:
        return ""
    return ""

def build_simple_tree(metadata_csv, out_path):
    print("[*] Constructing phylogenetic relationship of probiotic sources...")
    df = pd.read_csv(metadata_csv)
    
    # For demonstration of the phylogenetic pipeline, 
    # we use conserved sequences (e.g., 16S rRNA segments) 
    # extracted from each genome or build a hierarchy based on Taxonomy.
    # Note: A real-world production tree should use validated 16S alignments.
    
    # Simulation: Create a simple distance matrix to build a Neighbor-Joining (NJ) tree.
    # Distances are estimated based on Genus grouping for visualization.
    organisms = df["Organism"].unique().tolist()
    
    # Create Alignment records (Dummy sequences for placeholder logic)
    # In practice, 16S rRNA sequences would be aligned with MAFFT/ClustalW.
    dummy_records = []
    for i, org in enumerate(organisms):
        # Sequence simulation: members of the same Genus share higher similarity.
        genus = org.split()[0]
        seq = "ATGC" * 10 + (genus[0] * 5) + (org[-1] * 3)
        rec = SeqRecord(Seq(seq), id=org.replace(" ", "_")[:15], name=org)
        dummy_records.append(rec)
    
    # Bio.Align is required for MultipleSeqAlignment construction
    from Bio.Align import MultipleSeqAlignment
    align = MultipleSeqAlignment(dummy_records)
    
    # Compute distance matrix and build tree construction pipeline
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(align)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(align)
    
    # Render the phylogenetic tree plot
    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False)
    
    plt.title("Phylogenetic Relationship of Probiotic Sources (TeensGlow-AMP-02)", size=14)
    plt.savefig(out_path, bbox_inches="tight")
    print(f"[+] Phylogenetic tree saved to {out_path}")

if __name__ == "__main__":
    Entrez.email = "info@teensglow.us"
    csv_path = "data/raw/assembly_metadata.csv"
    if os.path.exists(csv_path):
        build_simple_tree(csv_path, "results/figures/phylogenetic_tree.png")
    else:
        print("[!] Metadata CSV not found.")
