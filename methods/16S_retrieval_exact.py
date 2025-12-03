#!/usr/bin/env python3

import argparse
import csv
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from urllib.error import HTTPError

# Global config
NCBI_SLEEP_BETWEEN_REQUESTS = 0.5
NCBI_MAX_RETRIES = 4
MIN_LEN = 100      # absolute min length to accept
MAX_LEN = 5000     # absolute max length to accept
AMBIGUOUS_MAX_FRAC = 0.5  # reject if >50% non-ACGT

taxonomy_cache: Dict[str, Dict] = {}


def configure_logging(log_file: Path):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)


def entrez_call(func, **kwargs):
    for attempt in range(1, NCBI_MAX_RETRIES + 1):
        try:
            handle = func(**kwargs)
            return handle
        except HTTPError as e:
            if e.code == 429 or e.code >= 500:
                wait = NCBI_SLEEP_BETWEEN_REQUESTS * (2 ** (attempt - 1))
                logging.warning(f"Entrez HTTP error {e.code}, retry {attempt}/{NCBI_MAX_RETRIES} after {wait:.1f}s")
                time.sleep(wait)
            else:
                logging.error(f"Entrez HTTP error {e.code}: {e}")
                raise
        except Exception as e:
            logging.warning(f"Entrez error on attempt {attempt}/{NCBI_MAX_RETRIES}: {e}")
            time.sleep(NCBI_SLEEP_BETWEEN_REQUESTS * (2 ** (attempt - 1)))
    raise RuntimeError("Entrez call failed after maximum retries")


def get_taxonomy_info(name: str) -> Dict:
    key = name.strip().lower()
    if key in taxonomy_cache:
        return taxonomy_cache[key]

    info = {
        "taxid": None,
        "canonical_name": name,
        "all_names": {name},
    }

    try:
        with entrez_call(Entrez.esearch, db="taxonomy", term=name, retmode="xml") as handle:
            data = Entrez.read(handle)
        idlist = data.get("IdList", [])
        if not idlist:
            logging.info(f"Taxonomy: no record found for '{name}'")
            taxonomy_cache[key] = info
            return info

        taxid = idlist[0]
        info["taxid"] = taxid

        with entrez_call(Entrez.efetch, db="taxonomy", id=taxid, retmode="xml") as handle:
            records = Entrez.read(handle)
        if not records:
            taxonomy_cache[key] = info
            return info

        rec = records[0]
        canon = rec.get("ScientificName", name)
        info["canonical_name"] = canon
        names: Set[str] = {name, canon}

        other = rec.get("OtherNames", {})
        for field in ("Synonym", "GenbankSynonym", "EquivalentName", "Name"):
            vals = other.get(field, [])
            if isinstance(vals, str):
                names.add(vals)
            else:
                for v in vals:
                    if isinstance(v, str):
                        names.add(v)
        info["all_names"] = {n.strip() for n in names if n.strip()}

        taxonomy_cache[key] = info
        return info
    except Exception as e:
        logging.warning(f"Taxonomy lookup failed for '{name}': {e}")
        taxonomy_cache[key] = info
        return info


def organism_matches_species(organism: str, target_species: str) -> bool:
    org = organism.strip().split()
    tgt = target_species.strip().split()
    if len(org) < 2 or len(tgt) < 2:
        return False

    genus_o, species_o = org[0].lower(), org[1].lower()
    genus_t, species_t = tgt[0].lower(), tgt[1].lower()

    if genus_o != genus_t:
        return False
    if species_o in {"sp", "sp.", "spp.", "spp"}:
        return False
    if species_o != species_t:
        return False

    return True


def is_16s_feature(feature) -> bool:
    if feature.type.lower() != "rRNA".lower() and feature.type.lower() != "rrna":
        return False

    prod = ""
    gene = ""
    if "product" in feature.qualifiers:
        prod = feature.qualifiers["product"][0].lower()
    if "gene" in feature.qualifiers:
        gene = feature.qualifiers["gene"][0].lower()

    text = prod + " " + gene

    if "12s" in text:
        return False

    if ("16s" in text) or ("16 s" in text):
        return True
    if "rrnl" in text:
        return True
    if "large subunit ribosomal rna" in text:
        return True

    return False


def is_sequence_plausible(seq: str) -> bool:
    s = seq.upper()
    length = len(s)
    if length < MIN_LEN or length > MAX_LEN:
        return False

    acgt = sum(1 for c in s if c in "ACGT")
    if acgt == 0:
        return False
    frac_acgt = acgt / length
    if frac_acgt < (1 - AMBIGUOUS_MAX_FRAC):
        return False

    return True


def score_candidate(seq: str, product: str, accession: str, from_mito: bool) -> float:
    length = len(seq)
    score = float(length)

    prod = (product or "").lower()

    if prod.strip() == "16s ribosomal rna":
        score += 200.0
    elif "16s" in prod:
        score += 100.0
    if "rrnl" in prod or "large subunit ribosomal rna" in prod:
        score += 80.0

    acc = accession.upper()
    if acc.startswith("NC_") or acc.startswith("NR_") or acc.startswith("NZ_"):
        score += 150.0

    if from_mito:
        score += 50.0

    acgt = sum(1 for c in seq.upper() if c in "ACGT")
    frac_acgt = acgt / max(len(seq), 1)
    score += frac_acgt * 50.0

    return score


def search_nuccore_ids(query: str, retmax: int = 10) -> List[str]:
    with entrez_call(Entrez.esearch, db="nuccore", term=query, retmax=retmax, retmode="xml") as handle:
        data = Entrez.read(handle)
    return data.get("IdList", [])


def fetch_genbank_records(ids: List[str]) -> List[SeqRecord]:
    if not ids:
        return []
    with entrez_call(Entrez.efetch, db="nuccore", id=",".join(ids), rettype="gb", retmode="text") as handle:
        records = list(SeqIO.parse(handle, "gb"))
    return records


def find_16s_for_species(scientific_name: str, name_variants: List[str]) -> Tuple[Optional[Dict], List[str]]:
    candidates: List[Dict] = []
    debug_lines: List[str] = []

    strategies = [
        ("MITOGENOME_FEATURE", '"{name}"[Organism] AND mitochondrion[filter] AND ("complete genome"[Title] OR "mitochondrial genome"[Title])'),
        ("DIRECT_16S_GENE", '"{name}"[Organism] AND (16S[Gene] OR "16S ribosomal RNA"[Title])'),
        ("RRNL_GENE", '"{name}"[Organism] AND (rrnL[Gene] OR "large subunit ribosomal RNA"[Title]) AND mitochondrion[filter]'),
        ("RIBOSOMAL_MITO", '"{name}"[Organism] AND ribosomal[Title] AND mitochondrion[filter]'),
    ]

    seen_accessions: Set[str] = set()

    for n in name_variants:
        debug_lines.append(f"  Trying name variant: {n}")
        for strat_name, q_template in strategies:
            query = q_template.format(name=n)
            debug_lines.append(f"    Strategy {strat_name}: {query}")
            try:
                ids = search_nuccore_ids(query, retmax=10)
            except Exception as e:
                msg = f"      NCBI search failed for strategy {strat_name} name '{n}': {e}"
                debug_lines.append(msg)
                logging.warning(msg)
                continue

            if not ids:
                debug_lines.append(f"      No nuccore hits")
                continue

            debug_lines.append(f"      Found {len(ids)} nuccore IDs")

            try:
                records = fetch_genbank_records(ids)
            except Exception as e:
                msg = f"      Failed to fetch GenBank records: {e}"
                debug_lines.append(msg)
                logging.warning(msg)
                continue

            for rec in records:
                acc = rec.id
                if acc in seen_accessions:
                    continue
                seen_accessions.add(acc)

                organism = rec.annotations.get("organism", "")
                if not organism_matches_species(organism, scientific_name):
                    debug_lines.append(f"        Skip accession {acc}: organism '{organism}' != '{scientific_name}'")
                    continue

                from_mito = False
                if "source" in [f.type for f in rec.features]:
                    for f in rec.features:
                        if f.type == "source":
                            if "organelle" in f.qualifiers:
                                if any("mitochondrion" in v.lower() for v in f.qualifiers["organelle"]):
                                    from_mito = True
                            break

                for feat in rec.features:
                    if not is_16s_feature(feat):
                        continue
                    loc = feat.location
                    seq = rec.seq[loc.start:loc.end]
                    if loc.strand == -1:
                        seq = seq.reverse_complement()
                    seq_str = str(seq)

                    prod = ""
                    gene = ""
                    if "product" in feat.qualifiers:
                        prod = feat.qualifiers["product"][0]
                    if "gene" in feat.qualifiers:
                        gene = feat.qualifiers["gene"][0]

                    if not is_sequence_plausible(seq_str):
                        debug_lines.append(
                            f"        Reject feature in {acc}: implausible length/sequence (len={len(seq_str)}, product='{prod}')"
                        )
                        continue

                    score = score_candidate(seq_str, prod, acc, from_mito)
                    candidate = {
                        "sequence": seq_str,
                        "accession": acc,
                        "product": prod,
                        "gene": gene,
                        "organism": organism,
                        "start": int(loc.start),
                        "end": int(loc.end),
                        "strand": int(loc.strand) if loc.strand is not None else 1,
                        "from_mito": from_mito,
                        "score": score,
                        "search_name": n,
                        "strategy": strat_name,
                    }
                    debug_lines.append(
                        f"        Accept candidate from {acc}: len={len(seq_str)}, score={score:.1f}, product='{prod}'"
                    )
                    candidates.append(candidate)

    if not candidates:
        return None, debug_lines

    candidates.sort(key=lambda c: c["score"], reverse=True)
    best = candidates[0]
    debug_lines.append(
        f"  BEST: accession={best['accession']}, len={len(best['sequence'])}, score={best['score']:.1f}, "
        f"product='{best['product']}', strategy={best['strategy']}, search_name='{best['search_name']}'"
    )

    return best, debug_lines


def detect_columns(header: List[str], scientific_column_arg: Optional[str]) -> Tuple[int, Optional[int]]:
    header_lower = [h.strip().lower() for h in header]

    sci_idx = None
    if scientific_column_arg:
        name_lower = scientific_column_arg.strip().lower()
        if name_lower in header_lower:
            sci_idx = header_lower.index(name_lower)
        else:
            raise ValueError(f"Scientific column '{scientific_column_arg}' not found in header: {header}")

    if sci_idx is None:
        candidates = ["scientific", "scientific_name", "sci_name", "species", "host_scientific_name"]
        for cand in candidates:
            for i, h in enumerate(header_lower):
                if cand == h or cand in h:
                    sci_idx = i
                    break
            if sci_idx is not None:
                break

    if sci_idx is None:
        if len(header) >= 2:
            sci_idx = 1
        else:
            sci_idx = 0

    common_idx = None
    for cand in ["common_name", "common", "host_common_name"]:
        if cand in header_lower:
            common_idx = header_lower.index(cand)
            break

    return sci_idx, common_idx


def load_species_list(path: Path, scientific_column: Optional[str]) -> List[Dict]:
    with path.open("r", encoding="utf-8", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.get_dialect("excel")

        reader = csv.reader(f, dialect)
        rows = list(reader)

    if not rows:
        raise ValueError("Input file appears to be empty")

    header = rows[0]
    sci_idx, common_idx = detect_columns(header, scientific_column)

    species_list: List[Dict] = []
    for idx, row in enumerate(rows[1:], start=1):
        if not row or all(not cell.strip() for cell in row):
            continue
        if sci_idx >= len(row):
            continue
        scientific_name = row[sci_idx].strip()
        if not scientific_name:
            continue
        common_name = row[common_idx].strip() if (common_idx is not None and common_idx < len(row)) else ""
        species_list.append(
            {
                "input_index": idx,
                "scientific_name": scientific_name,
                "common_name": common_name,
            }
        )

    return species_list


def main():
    parser = argparse.ArgumentParser(
        description="Exact 16S rRNA retriever: species-only, feature-based 16S/rrnL extraction from NCBI nuccore."
    )
    parser.add_argument("--email", required=True, help="Email address for NCBI Entrez")
    parser.add_argument("--api-key", default=None, help="NCBI API key (optional)")
    parser.add_argument("--input", required=True, help="Input CSV/TSV with at least scientific name column")
    parser.add_argument("--scientific-column", default=None, help="Header name of scientific name column (optional)")
    parser.add_argument("--output-fasta", default="combined_16S.fasta", help="Output FASTA file")
    parser.add_argument("--sources-csv", default="16S_sources.csv", help="Metadata CSV for successful sequences")
    parser.add_argument("--unfound-csv", default="unfound_species.csv", help="CSV listing species with no 16S found")
    parser.add_argument("--log-file", default="progress_log.txt", help="Log file with detailed progress")

    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    input_path = Path(args.input)
    fasta_path = Path(args.output_fasta)
    sources_path = Path(args.sources_csv)
    unfound_path = Path(args.unfound_csv)
    log_path = Path(args.log_file)

    configure_logging(log_path)

    logging.info(f"Input file: {input_path}")
    logging.info(f"Output FASTA: {fasta_path}")
    logging.info(f"Sources CSV: {sources_path}")
    logging.info(f"Unfound CSV: {unfound_path}")
    logging.info("Loading species list...")

    try:
        species_list = load_species_list(input_path, args.scientific_column)
    except Exception as e:
        logging.error(f"Failed to load species list: {e}")
        sys.exit(1)

    logging.info(f"Loaded {len(species_list)} species rows")

    fasta_out = fasta_path.open("w", encoding="utf-8")
    sources_out = sources_path.open("w", encoding="utf-8", newline="")
    unfound_out = unfound_path.open("w", encoding="utf-8", newline="")

    sources_writer = csv.writer(sources_out)
    sources_writer.writerow(
        [
            "input_index",
            "scientific_name",
            "common_name",
            "taxid",
            "final_name_used",
            "accession",
            "feature_start",
            "feature_end",
            "strand",
            "length",
            "database",
            "nuccore_url",
            "genbank_url",
            "search_strategy",
            "product",
            "gene",
        ]
    )

    unfound_writer = csv.writer(unfound_out)
    unfound_writer.writerow(
        [
            "input_index",
            "scientific_name",
            "common_name",
            "taxid",
            "names_tried",
            "reason",
            "notes",
        ]
    )

    for i, row in enumerate(species_list, start=1):
        input_index = row["input_index"]
        scientific_name = row["scientific_name"]
        common_name = row["common_name"]

        logging.info(f"[{i}/{len(species_list)}] Processing '{scientific_name}' (row {input_index})")

        tax_info = get_taxonomy_info(scientific_name)
        taxid = tax_info["taxid"]
        canonical_name = tax_info["canonical_name"]
        names_set = tax_info["all_names"]

        names_set.add(scientific_name)
        name_variants = sorted({n for n in names_set if n.strip()})

        logging.info(
            f"  TaxID: {taxid if taxid else 'N/A'}, canonical='{canonical_name}', "
            f"name_variants={'; '.join(name_variants)}"
        )

        best, debug_lines = find_16s_for_species(canonical_name, name_variants)

        for line in debug_lines:
            logging.info(line)

        if best is None:
            reason = "NO_16S_RRNA_FEATURE"
            notes = "No acceptable 16S/rrnL features found in any nuccore record"
            logging.warning(f"  FAILED for '{scientific_name}': {reason}")
            unfound_writer.writerow(
                [
                    input_index,
                    scientific_name,
                    common_name,
                    taxid if taxid else "",
                    "; ".join(name_variants),
                    reason,
                    notes,
                ]
            )
            unfound_out.flush()
            continue

        seq = best["sequence"]
        accession = best["accession"]
        start = best["start"]
        end = best["end"]
        strand = best["strand"]
        product = best["product"]
        gene = best["gene"]
        search_name = best["search_name"]
        strategy = best["strategy"]

        fasta_header = scientific_name
        fasta_out.write(f">{fasta_header}\n")
        for pos in range(0, len(seq), 60):
            fasta_out.write(seq[pos : pos + 60] + "\n")
        fasta_out.flush()

        db = "NCBI_nuccore"
        nuccore_url = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
        genbank_url = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}?report=genbank"
        length = len(seq)

        sources_writer.writerow(
            [
                input_index,
                scientific_name,
                common_name,
                taxid if taxid else "",
                search_name,
                accession,
                start,
                end,
                strand,
                length,
                db,
                nuccore_url,
                genbank_url,
                strategy,
                product,
                gene,
            ]
        )
        sources_out.flush()

        logging.info(
            f"  SUCCESS for '{scientific_name}': accession={accession}, len={length}, "
            f"product='{product}', strategy={strategy}, search_name='{search_name}'"
        )

    fasta_out.close()
    sources_out.close()
    unfound_out.close()
    logging.info("Done.")


if __name__ == "__main__":
    main()
