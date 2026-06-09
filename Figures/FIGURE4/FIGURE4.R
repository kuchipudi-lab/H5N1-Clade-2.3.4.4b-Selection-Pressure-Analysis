library(circlize)
library(dplyr)
library(readxl)

# ── DATA ──────────────────────────────────────────────────────────────────────

df <- read_excel(
  "C:/Users/Supplementary Table 1 (Updated).xlsx"
) %>%
  filter(!is.na(Gene))

GENE_ORDER <- c("PB2","PB1","PA","HA","NP","M1","M2","NS1","NS2")

GENE_LEN <- c(
  PB2 = 759,
  PB1 = 757,
  PA  = 716,
  HA  = 566,
  NP  = 498,
  M1  = 252,
  M2  = 97,
  NS1 = 230,
  NS2 = 121
)

GENE_COL <- c(
  PB2 = "#2B5C8D",
  PB1 = "#4183C4",
  PA  = "#74B2DB",
  HA  = "#C87535",
  NP  = "#D4B275",
  M1  = "#7450A0",
  M2  = "#9B70C8",
  NS1 = "#2F6848",
  NS2 = "#4EAA6A"
)

MC <- c(
  MEME  = "#1F77B4",
  cFEL  = "#FF7F0E",
  FUBAR = "#D62728",
  SLAC  = "#2CA02C",
  FEL   = "#9467BD"
)

# ── FLAGS ─────────────────────────────────────────────────────────────────────

df <- df %>%
  mutate(
    sig_MEME  = !is.na(MEME)               & MEME < 0.05,
    sig_cFEL  = !is.na(`Contrast-FEL`)     & `Contrast-FEL` < 0.05,
    sig_FUBAR = !is.na(FUBAR_post)         & FUBAR_post > 0.90,
    sig_SLAC  = !is.na(SLAC_p)             & SLAC_p < 0.05,
    sig_FEL   = !is.na(FEL_positive_p_val) & FEL_positive_p_val < 0.05,
    n_methods = sig_MEME + sig_cFEL + sig_FUBAR + sig_SLAC + sig_FEL,
    bold_label = sig_cFEL & n_methods >= 2
  ) %>%
  filter(Gene %in% GENE_ORDER) %>%
  mutate(Gene = factor(Gene, levels = GENE_ORDER))

cat("Total rows after Gene filter:", nrow(df), "\n")
cat("Sites labeled (>=2 methods):", sum(df$n_methods >= 2), "\n")
cat("Bold labels (Contrast-FEL + >=1 other):", sum(df$bold_label), "\n")

# ── MUTATION LABELS ───────────────────────────────────────────────────────────

parse_dom_aa <- function(s) {
  if (is.na(s) || s == "-") return("?")
  first <- trimws(strsplit(s, ";")[[1]][1])
  m <- regmatches(first, regexpr("^[A-Z*]", first))
  if (length(m) == 0) return("?")
  return(m)
}

df$dom_aa    <- sapply(df$Composition_2.3.4.4b, parse_dom_aa)
df$mut_label <- paste0(df$dom_aa, df$Codon)

marker_data <- df %>%
  filter(n_methods >= 2) %>%
  mutate(Gene = as.character(Gene)) %>%
  arrange(Gene, Codon)

bed <- marker_data %>%
  transmute(
    chr        = Gene,
    start      = as.numeric(Codon) - 0.5,
    end        = as.numeric(Codon) + 0.5,
    label      = mut_label,
    bold_label = bold_label
  ) %>%
  filter(chr %in% GENE_ORDER, start >= 0) %>%
  arrange(chr, start)

label_col  <- GENE_COL[bed$chr]
label_font <- ifelse(bed$bold_label, 2, 1)

# ── FONT SIZE CEX VALUES ──────────────────────────────────────────────────────
# Base font size in R graphics = 12 pt
# cex = target_pt / 12

CEX_GENE     <- 12 / 12   # 1.000 → 12 pt (max)
CEX_MUTATION <- 10 / 12   # 0.833 → 10 pt
CEX_TICK     <-  6 / 12   # 0.500 →  6 pt (min)
CEX_LEGEND   <- 10 / 12   # 0.833 → 10 pt

# ── TICK RING HELPER ──────────────────────────────────────────────────────────

draw_method_ticks <- function(gene, sig_col, col) {
  hits <- df$Codon[df$Gene == gene & df[[sig_col]]]
  for (p in hits) {
    circos.lines(
      c(p, p),
      c(0.08, 0.92),
      col = col,
      lwd = 2.0,
      straight = TRUE
    )
  }
}

# ── DRAWING FUNCTION ──────────────────────────────────────────────────────────

draw_circos_plot <- function() {
  
  circos.clear()
  par(mar = c(0.05, 0.05, 0.05, 0.05))
  
  circos.par(
    start.degree = 90,
    gap.after    = rep(2.5, length(GENE_ORDER)),
    clock.wise   = TRUE,
    # Tighter canvas so circle fills more of the 7x7 space
    # and labels still have room outside
    canvas.xlim  = c(-1.05, 1.05),
    canvas.ylim  = c(-1.05, 1.05),
    track.margin = c(0.004, 0.004),
    cell.padding = c(0, 0, 0, 0)
  )
  
  circos.initialize(
    factors = GENE_ORDER,
    xlim = cbind(
      rep(0, length(GENE_ORDER)),
      GENE_LEN[GENE_ORDER]
    )
  )
  
  # Mutation labels outside the circle — 10 pt
  circos.genomicLabels(
    bed,
    labels.column     = 4,
    side              = "outside",
    col               = rep("black", nrow(bed)),
    line_col          = label_col,
    cex               = CEX_MUTATION,
    font              = label_font,
    connection_height = mm_h(2.5),
    labels_height     = 0.18
  )
  
  # Gene arc — 12 pt
  circos.track(
    ylim         = c(0, 1),
    track.height = 0.14,
    bg.col       = NA,
    bg.border    = NA,
    panel.fun    = function(x, y) {
      
      gene <- get.cell.meta.data("sector.index")
      xl   <- get.cell.meta.data("xlim")
      plen <- GENE_LEN[gene]
      
      circos.rect(
        xl[1], 0, xl[2], 1,
        col    = GENE_COL[gene],
        border = "white",
        lwd    = 1.2
      )
      
      circos.text(
        plen / 2, 0.25, gene,
        facing = "bending.inside",
        font   = 2,
        cex    = CEX_GENE,
        col    = "white"
      )
      
      if (plen > 25) {
        
        # Major ticks every 100 aa — 6 pt
        if (plen > 100) {
          for (sp in seq(100, plen - 1, by = 100)) {
            circos.lines(
              c(sp, sp), c(0.75, 1.0),
              col = "white", lwd = 1.1, straight = TRUE
            )
            circos.text(
              sp, 0.63, as.character(sp),
              facing = "inside", niceFacing = TRUE,
              cex = CEX_TICK, col = "white", adj = c(0.5, 0.5)
            )
          }
        }
        
        # Mid ticks every 50 aa
        for (sp in seq(50, plen - 1, by = 50)) {
          if (sp %% 100 != 0) {
            circos.lines(
              c(sp, sp), c(0.84, 1.0),
              col = "white", lwd = 0.7, straight = TRUE
            )
          }
        }
        
        # Minor ticks every 25 aa
        for (sp in seq(25, plen - 1, by = 25)) {
          if (sp %% 50 != 0) {
            circos.lines(
              c(sp, sp), c(0.91, 1.0),
              col = "white", lwd = 0.5, straight = TRUE
            )
          }
        }
      }
    }
  )
  
  # Five method tick rings
  circos.track(
    ylim = c(0,1), track.height = 0.042,
    bg.col = "#f2f0ec", bg.border = NA,
    panel.fun = function(x,y) {
      draw_method_ticks(get.cell.meta.data("sector.index"), "sig_MEME", MC["MEME"])
    }
  )
  
  circos.track(
    ylim = c(0,1), track.height = 0.042,
    bg.col = "#f2f0ec", bg.border = NA,
    panel.fun = function(x,y) {
      draw_method_ticks(get.cell.meta.data("sector.index"), "sig_cFEL", MC["cFEL"])
    }
  )
  
  circos.track(
    ylim = c(0,1), track.height = 0.042,
    bg.col = "#f2f0ec", bg.border = NA,
    panel.fun = function(x,y) {
      draw_method_ticks(get.cell.meta.data("sector.index"), "sig_FUBAR", MC["FUBAR"])
    }
  )
  
  circos.track(
    ylim = c(0,1), track.height = 0.042,
    bg.col = "#f2f0ec", bg.border = NA,
    panel.fun = function(x,y) {
      draw_method_ticks(get.cell.meta.data("sector.index"), "sig_SLAC", MC["SLAC"])
    }
  )
  
  circos.track(
    ylim = c(0,1), track.height = 0.042,
    bg.col = "#f2f0ec", bg.border = NA,
    panel.fun = function(x,y) {
      draw_method_ticks(get.cell.meta.data("sector.index"), "sig_FEL", MC["FEL"])
    }
  )
  
  # Legend - placed inside the circle centre using canvas coordinates
  # canvas = c(-1.05, 1.05); legend sits in hollow centre
  legend(
    x = 0.55, y = 1.12,
    legend = c("MEME", "Contrast-FEL", "FUBAR", "SLAC", "FEL"),
    col    = unname(MC),
    lty    = 1,
    lwd    = 3,
    bty    = "o",
    cex    = CEX_LEGEND,
    xjust  = 0,
    yjust  = 1,
    xpd    = TRUE,
    title      = expression(italic("Tick rings (outer \u2192 inner)")),
    title.adj  = 0
  )
}

# ── SAVE PDF (7 × 7 in) ───────────────────────────────────────────────────────

pdf("FIGURE4.pdf", width = 7, height = 7)
draw_circos_plot()
dev.off()

# ── SAVE PNG (7 in × 600 dpi = 4200 px) ──────────────────────────────────────

png(
  "FIGURE4.png",
  width  = 4200,
  height = 4200,
  res    = 600
)
draw_circos_plot()
dev.off()

circos.clear()

message("Done: FIGURE4.pdf")
message("Done: FIGURE4.png")