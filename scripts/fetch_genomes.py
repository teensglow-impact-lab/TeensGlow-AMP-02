"""
TeensGlow-AMP-02 | Module 1: Genome Acquisition
================================================
Fetches complete probiotic genomes from NCBI Assembly via Entrez API.
Stores metadata in a local CSV registry for reproducibility.

Usage:
    python fetch_genomes.py --query "Lactobacillus plantarum" --limit 2 --outdir data/raw
"""

import os
import argparse
import pandas as pd
from Bio import Entrez
from tqdm import tqdm
import requests
import time

Entrez.email = "your_email@example.com"  # Required by NCBI policy

def search_assemblies(query, retmax=10):
    """Search NCBI Assembly for complete genome records matching the query."""
    print(f"[1/3] Searching NCBI Assembly: '{query}' (complete genomes only)...")
    handle = Entrez.esearch(
        db="assembly",
        term=f"{query} AND \"complete genome\"[filter]",
        retmax=retmax
    )
    record = Entrez.read(handle)
    handle.close()
    ids = record["IdList"]
    print(f"      Found {len(ids)} assembly record(s).")
    return ids

def get_assembly_details(assembly_ids):
    """Retrieve FTP paths and organism metadata for each Assembly ID."""
    print(f"[2/3] Fetching metadata for {len(assembly_ids)} assemblies...")
    details = []
    for aid in tqdm(assembly_ids, desc="Fetching metadata"):
        handle = Entrez.esummary(db="assembly", id=aid)
        record = Entrez.read(handle)
        handle.close()

        summary = record["DocumentSummarySet"]["DocumentSummary"][0]
        ftp_path = summary["FtpPath_GenBank"]

        if ftp_path:
            details.append({
                "ID": aid,
                "Organism": summary["Organism"],
                "Name": summary["AssemblyName"],
                "FTP": ftp_path
            })
        time.sleep(0.4)  # Respect NCBI rate limit (max 3 req/sec without API key)
    return details

def download_genome(ftp_path, target_dir):
    """Download genomic FASTA (.fna.gz) from NCBI FTP via HTTPS."""
    base_url = ftp_path.replace("ftp://", "https://")
    file_name = ftp_path.split("/")[-1] + "_genomic.fna.gz"
    download_url = f"{base_url}/{file_name}"
    local_path = os.path.join(target_dir, file_name)

    if os.path.exists(local_path):
        return local_path  # Skip if already downloaded

    try:
        response = requests.get(download_url, stream=True, timeout=60)
        response.raise_for_status()
        with open(local_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return local_path
    except Exception as e:
        print(f"  [WARNING] Download failed for {download_url}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="TeensGlow-AMP-02: Probiotic Genome Downloader"
    )
    parser.add_argument("--query", type=str, required=True,
                        help="NCBI search query (e.g. 'Lactobacillus plantarum')")
    parser.add_argument("--limit", type=int, default=5,
                        help="Maximum number of genomes to download")
    parser.add_argument("--outdir", type=str, default="data/raw",
                        help="Output directory for downloaded genomes")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ids = search_assemblies(args.query, retmax=args.limit)
    if not ids:
        print("[ERROR] No assemblies found. Adjust query or check NCBI connectivity.")
        return

    details = get_assembly_details(ids)

    print(f"[3/3] Downloading {len(details)} genome(s)...")
    results = []
    for item in tqdm(details, desc="Downloading genomes"):
        path = download_genome(item["FTP"], args.outdir)
        if path:
            item["LocalPath"] = path
            results.append(item)

    # Update metadata registry (append, deduplicate by Assembly ID)
    csv_path = os.path.join(args.outdir, "assembly_metadata.csv")
    df_new = pd.DataFrame(results)
    if os.path.exists(csv_path):
        df_old = pd.read_csv(csv_path)
        df_new = pd.concat([df_old, df_new]).drop_duplicates(subset=["ID"])
    df_new.to_csv(csv_path, index=False)
    print(f"[DONE] Metadata registry saved to: {csv_path}")

if __name__ == "__main__":
    main()
