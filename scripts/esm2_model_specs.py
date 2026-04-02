"""
TeensGlow-AMP-02: ESM-2 Feature Extraction Details
Provides detailed model configuration and embedding dimensionalities 
as required for the Methods section of Computational Biology and Chemistry.
"""

import torch
from transformers import AutoTokenizer, EsmModel
import pandas as pd
import numpy as np

def get_model_specs():
    """
    Prints detailed ESM-2 model parameters for inclusion in the Methods section.
    """
    model_name = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)

    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

    specs = {
        "Model": model_name,
        "Architecture": "ESM-2 (Evolutionary Scale Modeling)",
        "Transformer Layers": 6,
        "Attention Heads": 20,
        "Embedding Dim": 320,
        "Total Parameters": f"{total_params:,}",
        "Training Corpus": "UR50/D (UniRef50 + UniRef50/D), 65M protein sequences",
        "Usage": "Zero-shot feature extraction (no fine-tuning)"
    }

    print("=" * 60)
    print("ESM-2 Model Specifications for Methods Section")
    print("=" * 60)
    for k, v in specs.items():
        print(f"  {k:25s}: {v}")
    print("=" * 60)

    # Construct example embedding to demonstrate dimensionality
    test_seq = "MRLLKLVIHLVKKLRLLLTRRHR"
    inputs = tokenizer(test_seq, return_tensors="pt", add_special_tokens=False)
    with torch.no_grad():
        output = model(**inputs)
    emb = output.last_hidden_state

    print(f"\n  Input sequence length : {len(test_seq)} AA")
    print(f"  Embedding shape       : {list(emb.shape)} (batch, seq_len, dim)")
    print(f"  Mean-pooled vector    : shape {list(emb.mean(dim=1).shape)}")

    # Save detailed parameters to CSV for supplementary materials
    df = pd.DataFrame([specs])
    df.to_csv("results/esm2_model_specs.csv", index=False)
    print(f"\n[+] Model specs saved to results/esm2_model_specs.csv (for Supplementary Materials)")
    return specs

if __name__ == "__main__":
    get_model_specs()
