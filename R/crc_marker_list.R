# Canonical CRC / TME marker gene list used for cluster annotation dot plots.
# Sourced by P*_seurat/S1.4_*_select_epi.R and
# seurat_integration/S2.2_intigre_after_selected.R (run from repo root).

crc_marker_list <- list(

  # 1. CRC tumor epithelial
  "CRC tumor epithelial" = c(
    "CEACAM5",
    "CEACAM6",
    "MUC17"
  ),

  # 2. Stem / progenitor epithelial
  "Stem-like epithelial" = c(
    "LGR5",
    "ASCL2",
    "OLFM4"
  ),

  # 3. BEST4+ / absorptive epithelial
  "BEST4+ epithelial" = c(
    "BEST4",
    "GUCA2A",
    "OTOP2"
  ),

  # 4. Secretory / goblet-like epithelial
  "Secretory epithelial" = c(
    "REG4",
    "REG1A",
    "CLCA1"
  ),

  # 5. Proliferating cells
  "Proliferating" = c(
    "MKI67",
    "UBE2C",
    "STMN1"
  ),

  # 6. Macrophage
  "Macrophage" = c(
    "C1QA",
    "C1QB",
    "C1QC"
  ),

  # 7. Monocyte / inflammatory myeloid
  "Monocyte inflammatory" = c(
    "FCN1",
    "S100A12",
    "VCAN"
  ),

  # 8. T cell
  "T cell" = c(
    "CD3D",
    "CD3E",
    "TRAC"
  ),

  # 9. Cytotoxic T / NK
  "Cytotoxic T NK" = c(
    "NKG7",
    "GNLY",
    "PRF1"
  ),

  # 10. Treg / activated T
  "Treg activated T" = c(
    "FOXP3",
    "IL2RA",
    "CTLA4"
  ),

  # 11. B cell
  "B cell" = c(
    "MS4A1",
    "CD79A",
    "PAX5"
  ),

  # 12. Plasma cell
  "Plasma cell" = c(
    "MZB1",
    "TNFRSF17",
    "FKBP11"
  ),

  # 13. CAF
  "CAF" = c(
    "PECAM1",
    "VWF",
    "CDH5"
  ),

  # 14. CAF / stromal
  "CAF stromal" = c(
    "FAP",
    "COL11A1",
    "CTHRC1"
  ),

  # 15. Pericyte / smooth muscle
  "Pericyte smooth muscle" = c(
    "RGS5",
    "ACTA2",
    "TAGLN"
  )
)
