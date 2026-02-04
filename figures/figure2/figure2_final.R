###############################################
## Fig 1 – Influenza mutation circos (PDF)   ##
## MEME + Contrast-FEL + BOTH (bold labels)  ##
###############################################

## 1. Libraries ----
library(circlize)
library(readxl)
library(dplyr)
library(tidyr)
library(grid)     # for central legend

## 2. Segment order, lengths, colours ----

# segment order around the circle
seg_order <- c("PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NS2")

# segment lengths (AA) in same order as seg_order
seg_lengths_values <- c(759, 757, 716, 566, 498, 469, 252, 97, 230, 121)
seg_lengths <- setNames(seg_lengths_values, seg_order)

# colour palette: names + values
seg_names  <- c("PB2","PB1","PA","HA","NA","NP","M1","M2","NS1","NS2")
seg_colors <- c(
  "#003049",  # PB2 – indigo
  "#264653",  # PB1 – slate
  "#006D5B",  # PA  – forest
  "#BB3E03",  # HA  – brick
  "#BC6C25",  # NA  – clay
  "#8FB339",  # NP  – moss
  "#6A4C93",  # M1  – plum
  "#6A4C93",  # M2  – same as M1
  "#006D5B",  # NS1 – forest
  "#006D5B"   # NS2 – forest
)
line_colors <- setNames(seg_colors, seg_names)
seg_cols    <- line_colors   # use same palette for ring fill

## 3. Read combined MEME + Contrast-FEL XLSX ----
results_path <- "MEME and Contrast FEL Results (Cleaned) .xlsx"

# pattern: AA + number (+ optional AA), e.g., A29, N71S
mutation_pattern <- "^[A-Z][0-9]+[A-Z]?$"

# read raw
df_raw <- read_excel(results_path)

# quick sanity check
if (!all(c("Segment", "MEME Position", "CONTRAST FEL Position", "BOTH") %in% names(df_raw))) {
  stop("Expected columns 'Segment', 'MEME Position', 'CONTRAST FEL Position', 'BOTH' not all found.")
}

## 3a. Long format: one row per (Segment, mutation, category) ----
mut_long_raw <- df_raw %>%
  pivot_longer(
    cols      = c(`MEME Position`, `CONTRAST FEL Position`, BOTH),
    names_to  = "Category_raw",
    values_to = "mutation"
  ) %>%
  mutate(
    Segment  = toupper(trimws(Segment)),
    mutation = trimws(as.character(mutation))
  ) %>%
  # keep only non-empty, valid AA+number(+AA) strings
  filter(!is.na(mutation),
         mutation != "",
         grepl(mutation_pattern, mutation)) %>%
  mutate(
    cat_flag = case_when(
      Category_raw == "MEME Position"         ~ "MEME",
      Category_raw == "CONTRAST FEL Position" ~ "CONTRAST-FEL",
      Category_raw == "BOTH"                  ~ "BOTH",
      TRUE                                    ~ "MEME"
    )
  )

## 3b. Collapse duplicates per site, prioritising BOTH ----
mut_long <- mut_long_raw %>%
  group_by(Segment, mutation) %>%
  summarise(
    has_meme = any(cat_flag == "MEME"),
    has_cfel = any(cat_flag == "CONTRAST-FEL"),
    has_both = any(cat_flag == "BOTH"),
    .groups  = "drop"
  ) %>%
  mutate(
    category = case_when(
      has_both ~ "BOTH",
      has_meme & has_cfel ~ "BOTH",      # if present in both even without explicit BOTH column
      has_meme ~ "MEME",
      has_cfel ~ "CONTRAST-FEL",
      TRUE     ~ "MEME"
    ),
    # numeric position from the mutation string
    pos_raw = as.numeric(gsub("[^0-9]", "", mutation))
  ) %>%
  # join segment lengths so we can clip to valid range
  left_join(
    data.frame(
      segment = seg_order,
      seg_len = as.numeric(seg_lengths),
      stringsAsFactors = FALSE
    ),
    by = c("Segment" = "segment")
  ) %>%
  mutate(
    pos     = pmin(pos_raw, seg_len),
    start   = pmax(0, pos - 0.5),
    end     = pmin(seg_len, pos + 0.5),
    segment = factor(Segment, levels = seg_order),
    category = factor(category,
                      levels = c("MEME", "CONTRAST-FEL", "BOTH"))
  )

# bed data for circlize: sector, start, end, label, category
bed <- mut_long %>%
  dplyr::select(
    sector   = segment,
    start,
    end,
    label    = mutation,
    category
  )

## 4. Open vector PDF device ----
pdf("Fig1_influenza_circos_MEME_CFEL.pdf",
    width = 7.1, height = 7.1, family = "sans")

## 5. Circos setup ----
circos.clear()

circos.par(
  start.degree   = 90,
  gap.degree     = rep(2, length(seg_order)),
  track.margin   = c(0.01, 0.01),
  cell.padding   = c(0, 0, 0, 0),
  canvas.xlim    = c(-1.1, 1.1),
  canvas.ylim    = c(-1.1, 1.1)
)

circos.initialize(
  factors = seg_order,
  xlim    = cbind(rep(0, length(seg_order)), seg_lengths[seg_order])
)

## 6. COLORED SEGMENT RING + SCALE JUST OUTSIDE ----
circos.trackPlotRegion(
  ylim         = c(0, 1),
  track.height = 0.10,
  bg.border    = NA,
  panel.fun    = function(x, y) {
    
    seg <- CELL_META$sector.index
    len <- seg_lengths[seg]
    
    # coloured ring band
    circos.rect(
      CELL_META$xlim[1], 0,
      CELL_META$xlim[2], 1,
      col    = seg_cols[seg],
      border = "white"
    )
    
    # segment label
    circos.text(
      mean(CELL_META$xlim), 0.5, seg,
      cex = 0.7, col = "white", font = 2
    )
    
    # --- TICKS & SCALE (special case NS1 to avoid overlap) ---
    at <- c(50, 200, 400, 600, len)
    at <- at[at <= len]
    
    if (seg == "NS1") {
      at      <- c(len)    # just one tick at the end
      labels  <- at
      lab_cex <- 0.30
    } else {
      labels  <- at
      lab_cex <- 0.35
    }
    
    circos.axis(
      h                 = "top",
      major.at          = at,
      labels            = labels,
      labels.cex        = lab_cex,
      labels.niceFacing = TRUE,
      lwd               = 0.4,
      col               = "#444444"
    )
  }
)

## 7. LABELS INSIDE – single track, style by category ----
##   - MEME-only: solid connectors, normal font
##   - Contrast-FEL-only: dashed connectors, normal font
##   - BOTH: solid connectors, **bold labels**

line_lty_vec <- ifelse(bed$category == "CONTRAST-FEL", 3, 1)   # dashed for CFEL
font_vec     <- ifelse(bed$category == "BOTH", 2, 1)           # bold for BOTH

circos.genomicLabels(
  bed,
  labels.column     = "label",
  side              = "inside",
  line_col          = line_colors[as.character(bed$sector)],
  line_lwd          = 0.7,
  line_lty          = line_lty_vec,
  col               = "black",
  font              = font_vec,
  cex               = 0.5,
  connection_height = mm_h(2),
  labels_height     = 0.55
)

## 8. CENTRAL LEGEND — TABLE-LIKE WITH COLUMNS ----

# Y positions (compact, evenly spaced)
title_y <- 0.60
y1 <- 0.52   # MEME row
y2 <- 0.46   # CFEL row
y3 <- 0.40   # BOTH row

# Column positions
line_x1 <- 0.36   # start of line (column 1)
line_x2 <- 0.44   # end of line (column 1)
text_x  <- 0.46   # start of text (column 2, right after line)

# Title centered
grid::grid.text(
  "Selection pressures",
  x = 0.50, y = title_y,
  gp = grid::gpar(fontsize = 10, fontface = "bold")
)

# Row 1: MEME – solid line + text in next "column"
grid::grid.lines(
  x = grid::unit(c(line_x1, line_x2), "npc"),
  y = grid::unit(c(y1, y1), "npc"),
  gp = grid::gpar(lwd = 1.1, lty = 1)
)
grid::grid.text(
  "episodic diversifying (MEME)",
  x = text_x, y = y1,
  just = "left",
  gp = grid::gpar(fontsize = 8)
)

# Row 2: Contrast-FEL – dashed line + text in same text column
grid::grid.lines(
  x = grid::unit(c(line_x1, line_x2), "npc"),
  y = grid::unit(c(y2, y2), "npc"),
  gp = grid::gpar(lwd = 1.1, lty = 3)
)
grid::grid.text(
  "intensified selection (Contrast-FEL)",
  x = text_x, y = y2,
  just = "left",
  gp = grid::gpar(fontsize = 8)
)

# Row 3: BOTH – bold line + bold text in same column
grid::grid.lines(
  x = grid::unit(c(line_x1, line_x2), "npc"),
  y = grid::unit(c(y3, y3), "npc"),
  gp = grid::gpar(lwd = 1.3, lty = 1)
)
grid::grid.text(
  "under both pressures (BOTH)",
  x = text_x, y = y3,
  just = "left",
  gp = grid::gpar(fontsize = 8, fontface = "bold")
)

## 9. Close device & clear ----
dev.off()
circos.clear()
