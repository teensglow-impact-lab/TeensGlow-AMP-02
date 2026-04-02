# TeensGlow-AMP-02: Targeted Anti-Acne Bacteriocin Discovery

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![AI Model: ESM-2](https://img.shields.io/badge/AI-ESM--2-brightgreen)](https://huggingface.co/facebook/esm2_t6_8M_UR50D)
[![Resolution: 300 DPI](https://img.shields.io/badge/Resolution-300_DPI-orange)](results/figures/)

## AI-Driven Discovery of Probiotic-Derived Anti-Acne Bacteriocins

**TeensGlow-AMP-02** is a high-fidelity computational project for identifying and optimizing targeted antimicrobial peptides (AMPs) from probiotic genomes. By leveraging **ESM-2 transformer embeddings** and a proprietary **GlowScore** vetting function, we demonstrate a specialized screening pipeline for *Cutibacterium acnes* inhibitors with high sebum compatibility and minimal hemolysis.

This repository is optimized for manuscript submission to **Computational Biology and Chemistry (Elsevier)**.

---

## 🚀 Key Features

1. **ESM-2 AI Screening**: Zero-shot feature extraction using Facebook's `esm2_t6_8M_UR50D` transformer model (6 layers, 320-dim embeddings) to rank 4,000+ bacteriocin-like candidates from NCBI genomes.
2. **GlowScore Algorithm**: A multi-property scoring function integrating Antibacterial Potential, Targeting Index (TI), and Sebum-Specific Compatibility for adolescent skin conditions.
3. **Advanced Truncation & Docking**: Sliding-window Alpha-optimization to refine 20-30 AA fragments, validated by AutoDock Vina molecular docking against the *C. acnes* PBP2 target.
4. **Publication Figure Suite**: 300 DPI figures (including Graphical Abstract) meeting Elsevier scientific reporting standards.

## 🧬 Lead Candidate: Alpha-Core

- **Sequence**: `MRLLKLVIHLVKKLRLLLTRRHR` (23 AA)
- **Source**: *Lactiplantibacillus plantarum* genomic ORF.
- **Activity Profile**: High Targeting Index (TI=5.68), low hemolysis (<25%), and superior predicted binding energy ($\Delta G = -8.08$ kcal/mol).
- **Target**: *Cutibacterium acnes* Penicillin-Binding Protein 2 (PBP2).

---

## 📂 Repository Structure

- `scripts/`: Implementation of the screening pipeline, docking simulation, and figure generation.
- `data/`: Contains probiotic genome metadata and processed ORF libraries.
- `results/`: 
  - `figures/`: High-resolution figures (Fig 1-4, Graphical Abstract) for manuscript integration.
  - `supplements/`: CSV results tables including model specifications and docking energies.
- `CHANGELOG.md`: Detailed development history and publication milestones.

## 🛠️ Getting Started

```bash
# 1. Install dependencies
pip install torch transformers pandas biopython matplotlib seaborn scikit-learn

# 2. Extract model specifications
python scripts/esm2_model_specs.py

# 3. Generate all publication figures
python scripts/generate_pub_figures.py
```

## 📝 SCI Citation / Journal Reference

This work is part of the **TeensGlow-AMP-02** platform. 

> *Title: AI-Driven Discovery of Probiotic-Derived Anti-Acne Bacteriocins via ESM-2 Screening and Molecular Docking Analysis.*

**Corresponding Author**: [info@teensglow.us](mailto:info@teensglow.us) (TeensGlow Impact Lab)  
**Target Journal**: *Computational Biology and Chemistry (CBC)*  
**Current Status**: Publication Ready (v1.0.0)

---
© 2026 TeensGlow-AMP-02 Project. All Rights Reserved.
