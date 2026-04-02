# Changelog: TeensGlow-AMP-02

All notable changes to the **TeensGlow-AMP-02** project will be documented in this file.
This repository is prepared for submission to *Computational Biology and Chemistry (Elsevier)*.

## [1.0.0] - 2026-04-02 (Current / Publication Ready)
### Added
- **Final Figure Suite**: Generated five 300-DPI, Elsevier-standard publication figures (Fig 1-4 and Graphical Abstract).
- **ESM-2 Model Specs**: Added `scripts/esm2_model_specs.py` to provide model architecture details for the Methods section.
- **Statistical Rigor**: Integrated Welch's t-test and significance markers into comparative plots.
- **Phylogenetic Analysis**: Added `scripts/phylogenetic_tree.py` for probiotic source alignment and clustering.
- **Enhanced Documentation**: Fully translated code comments, gitignore, and metadata to English for international submission.

### Changed
- **Graphical Abstract**: Redesigned from scratch with a 9x5 inch canvas to meet CBC mandatory requirements.
- **Workflow Optimization**: Refined the Alpha-Core sliding window truncation logic for better sebum compatibility.
- **Docking Simulation**: Updated Vina scoring parameters for *C. acnes* PBP2 specific targeting.

## [0.5.0] - 2026-03-28 (AI Model Upgrade)
### Added
- **ESM-2 Screening**: Replaced legacy peptide embeddings with `facebook/esm2_t6_8M_UR50D` for high-fidelity zero-shot screening.
- **GlowScore v2.1**: Integrated Targeting Index (TI) and Sebum Compatibility into the scoring function.

## [0.1.0] - 2026-03-24 (Project Kickoff)
### Added
- Initial genome mining pipeline for *L. plantarum*, *S. epidermidis*, and *C. acnes*.
- Basic ORF extraction and hemolysis prediction modules.

---
*Maintained by the TeensGlow Impact Lab ([info@teensglow.us](mailto:info@teensglow.us)).*
