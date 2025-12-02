###############################################
## Fig 1 – Influenza mutation circos (PDF)   ##
## Segment-colored ring + connectors         ##
###############################################

## 1. Libraries ----
library(circlize)
library(readxl)
library(dplyr)
library(tidyr)

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

## 3. Read MEME results & build mutation table ----

# all mutations are in MEME_Results.xlsx
meme_path <- "MEME_Results.xlsx"

# pattern: AA + number (+ optional AA), e.g., A29, N71S
mutation_pattern <- "^[A-Z][0-9]+[A-Z]?$"

mut_long <- read_excel(meme_path) %>%
  # normalise and clean
  mutate(
    Segment  = if_else(is.na(Segment) | Segment == "", "NA", Segment),
    Segment  = toupper(trimws(Segment)),
    mutation = trimws(Position)
  ) %>%
  # keep only real mutation strings like A29, N71S, etc.
  filter(!is.na(mutation),
         grepl(mutation_pattern, mutation)) %>%
  # parse numeric position
  mutate(
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
  # clip to valid [1, seg_len]
  mutate(
    pos     = pmin(pos_raw, seg_len),
    start   = pmax(0, pos - 0.5),
    end     = pmin(seg_len, pos + 0.5),
    segment = factor(Segment, levels = seg_order)
  )

# bed data for circlize: sector, start, end, label
bed <- mut_long %>%
  dplyr::select(sector = segment, start, end, label = mutation)

## 4. Open vector PDF device ----
# ~7.1 in ≈ 180 mm
pdf("Fig1_influenza_circos.pdf",
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
      # minimal, non-overlapping scale for NS1
      at      <- c(len)    # just one tick at the end
      labels  <- at        # or rep("", length(at)) if you want no numbers
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

## 7. LABELS INSIDE WITH SEGMENT-COLORED CONNECTORS ----
circos.genomicLabels(
  bed,
  labels.column     = "label",
  side              = "inside",
  line_col          = line_colors[as.character(bed$sector)],  # connectors in segment colours
  line_lwd          = 0.7,
  col               = "black",      # text colour
  cex               = 0.5,
  connection_height = mm_h(2),
  labels_height     = 0.55
)

## 8. Close device & clear ----
dev.off()
circos.clear()
