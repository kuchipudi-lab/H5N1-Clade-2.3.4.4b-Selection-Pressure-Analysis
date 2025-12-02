#!/usr/bin/env python3
"""
Rescue mitochondrial 16S rRNA sequences for missing species from NCBI (nuccore).

Inputs (same directory):
  - combined_16S.fasta         : existing 16S dataset (258 sequences, unique)
  - missing_120_species.txt    : one scientific name per line (120 missing spp)

Outputs:
  - rescue_16S_from_ncbi.fasta : newly rescued 16S sequences only
  - rescue_16S_log.csv         : detailed log for each attempted species
  - rescue_16S_no_hits.txt     : species for which no acceptable 16S was found
  - rescue_16S_duplicates_only.csv : species where only duplicate sequences
                                     (already present in combined_16S) were found

Search strategy:
  1. Try direct species-level search (including subspecies + 2-word species-level).
  2. If no valid 16S, try genus-level search for stand-ins (sister species).
  3. Stop when the first valid 16S rRNA feature is found.
  4. Reject any candidate whose nucleotide sequence is already present in:
        - combined_16S.fasta, OR
        - previously accepted rescued sequences.
"""

import time
import csv
from pathlib import Path
from typing import List, Optional, Tuple, Set, Dict

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord


# ------------------ CONFIG ------------------ #

ENTREZ_EMAIL = "yoj48@pitt.edu"  # <<< EDIT THIS
Entrez.email = ENTREZ_EMAIL
Entrez.tool = "rescue_16S_ncbi_script"

DELAY_SEC = 0.5  # delay between NCBI calls (seconds)

COMBINED_FASTA = "combined_16S.fasta"
MISSING_TXT = "missing_120_species.txt"

OUT_FASTA = "rescue_16S_from_ncbi.fasta"
OUT_LOG = "rescue_16S_log.csv"
OUT_NO_HITS = "rescue_16S_no_hits.txt"
OUT_DUP_ONLY = "rescue_16S_duplicates_only.csv"


# ------------------ HELPERS ------------------ #

def load_existing_sequences(
    fasta_path: str,
) -> Tuple[Set[str], Dict[str, Set[str]]]:
    """
    Load existing sequences and map seq_string -> set of species names.

    Returns:
        (seq_strings_set, seq_to_species_map)
    """
    seqs: Set[str] = set()
    seq_to_species: Dict[str, Set[str]] = {}

    if not Path(fasta_path).exists():
        print(f"[WARN] {fasta_path} not found; treating as empty.")
        return seqs, seq_to_species

    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq_str = str(rec.seq).upper().replace("-", "")
        desc = rec.description.strip()
        seqs.add(seq_str)
        seq_to_species.setdefault(seq_str, set()).add(desc)

    print(f"[INFO] Loaded {len(seqs)} existing sequences from {fasta_path}")
    return seqs, seq_to_species


def load_missing_species(txt_path: str) -> List[str]:
    """Load missing species names from text file."""
    names: List[str] = []
    with open(txt_path, "r") as f:
        for line in f:
            name = line.strip()
            if name:
                names.append(name)
    print(f"[INFO] Loaded {len(names)} missing species from {txt_path}")
    return names


def extract_16s_feature(record: SeqRecord) -> Optional[Tuple[str, str, str]]:
    """
    Look through rRNA features and extract a 16S / rrnL / large subunit region.

    Returns:
        (sequence_string, product, location_str) or None if not found.
    """
    for feat in record.features:
        if feat.type.lower() != "rrna":
            continue

        products = feat.qualifiers.get("product", [""])
        product_raw = products[0] if products else ""
        product = product_raw.lower()

        if (
            "16s" in product
            or "rrnl" in product
            or "large subunit ribosomal rna" in product
        ):
            try:
                sub_seq = feat.location.extract(record.seq)
                seq_str = str(sub_seq).upper()
                # Skip tiny junk fragments
                if len(seq_str) < 50:
                    continue
                loc_str = str(feat.location)
                return seq_str, product_raw, loc_str
            except Exception as e:
                print(f"[WARN] Failed to extract 16S feature: {e}")
                continue

    return None


def entrez_search_nuccore(term: str, retmax: int = 50) -> List[str]:
    """Search NCBI nuccore, return list of IDs."""
    handle = Entrez.esearch(db="nuccore", term=term, retmax=retmax)
    time.sleep(DELAY_SEC)
    result = Entrez.read(handle)
    handle.close()
    return result.get("IdList", [])


def entrez_fetch_gb(nuccore_id: str) -> Optional[SeqRecord]:
    """Fetch a GenBank record by nuccore ID."""
    try:
        handle = Entrez.efetch(
            db="nuccore", id=nuccore_id, rettype="gb", retmode="text"
        )
        time.sleep(DELAY_SEC)
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return record
    except Exception as e:
        print(f"[ERROR] efetch failed for {nuccore_id}: {e}")
        return None


def build_species_search_terms(scientific_name: str) -> List[str]:
    """
    Build a list of search terms for a given species name.

    Strategy:
      - full name
      - full name with mitochondrion filter
      - two-word species-level name (if subspecies),
        with and without mitochondrion filter
    """
    sci = scientific_name.strip()
    parts = sci.split()
    terms: List[str] = []

    if len(parts) >= 1:
        full = " ".join(parts)
        terms.append(f"\"{full}\"[Organism] AND mitochondrion[Filter]")
        terms.append(f"\"{full}\"[Organism]")

    if len(parts) >= 2:
        first_two = " ".join(parts[:2])
        if first_two != sci:
            terms.append(f"\"{first_two}\"[Organism] AND mitochondrion[Filter]")
            terms.append(f"\"{first_two}\"[Organism]")

    # De-duplicate while preserving order
    seen = set()
    unique_terms: List[str] = []
    for t in terms:
        if t not in seen:
            seen.add(t)
            unique_terms.append(t)

    return unique_terms


def build_genus_search_term(scientific_name: str) -> Optional[str]:
    """Return a genus-level search term for fallback."""
    parts = scientific_name.strip().split()
    if not parts:
        return None
    genus = parts[0]
    return f"\"{genus}\"[Organism] AND mitochondrion[Filter]"


def make_header(original: str, mapped: str) -> str:
    """Create FASTA header 'Original_from_Mapped' with underscores."""
    o = original.replace(" ", "_")
    m = mapped.replace(" ", "_")
    return f"{o}_from_{m}"


# ------------------ CORE SEARCH LOGIC ------------------ #

def try_direct_species(
    scientific_name: str,
    existing_seq_strings: Set[str],
    existing_seq_map: Dict[str, Set[str]],
) -> Tuple[Optional[SeqRecord], dict]:
    """
    Try to find a 16S rRNA for the exact species (or its 2-word species-level).

    Returns:
        (SeqRecord or None, log_dict)
    """
    log = {
        "original_species": scientific_name,
        "mapped_species": scientific_name,
        "accession": "",
        "feature_location": "",
        "product": "",
        "sequence_length": "",
        "search_level": "DIRECT",
        "status": "",
        "ncbi_url": "",
        "notes": "",
    }

    duplicates_seen: Set[str] = set()

    terms = build_species_search_terms(scientific_name)
    print(f"[INFO] DIRECT search terms for '{scientific_name}':")
    for t in terms:
        print(f"       - {t}")

    for term in terms:
        try:
            ids = entrez_search_nuccore(term, retmax=50)
        except Exception as e:
            log["status"] = "ERROR"
            log["notes"] = f"esearch failed: {e}"
            print(
                f"[ERROR] esearch failed for '{scientific_name}' "
                f"term '{term}': {e}"
            )
            continue

        if not ids:
            continue

        for nid in ids:
            record = entrez_fetch_gb(nid)
            if record is None:
                continue

            # Accept any record that has a proper 16S feature
            extracted = extract_16s_feature(record)
            if extracted is None:
                continue

            seq_str, product, loc_str = extracted
            if seq_str in existing_seq_strings:
                existing_species = existing_seq_map.get(seq_str, set())
                duplicates_seen.update(existing_species)
                print(
                    f"[INFO] DIRECT candidate for '{scientific_name}' "
                    f"is a duplicate of existing: {', '.join(existing_species)}; skipping."
                )
                continue

            accession = record.id
            header = make_header(scientific_name, scientific_name)
            new_rec = SeqRecord(
                record.seq.__class__(seq_str),
                id=header,
                name=header,
                description=header,
            )

            existing_seq_strings.add(seq_str)

            log["accession"] = accession
            log["feature_location"] = loc_str
            log["product"] = product
            log["sequence_length"] = str(len(seq_str))
            log["status"] = "DIRECT_16S_FOUND"
            log["ncbi_url"] = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
            log["notes"] = ""
            return new_rec, log

    if duplicates_seen:
        log["status"] = "ONLY_DUPLICATE_16S_DIRECT"
        log["notes"] = (
            "Only duplicate 16S sequences found at species level; existing species: "
            + ", ".join(sorted(duplicates_seen))
        )
    else:
        log["status"] = "NO_DIRECT_16S"
        log["notes"] = (
            "No acceptable 16S found at species level "
            "(no 16S features or no records)."
        )

    return None, log


def try_genus_standin(
    scientific_name: str,
    existing_seq_strings: Set[str],
    existing_seq_map: Dict[str, Set[str]],
) -> Tuple[Optional[SeqRecord], dict]:
    """
    Try to find a 16S rRNA in the same genus as a stand-in.

    Returns:
        (SeqRecord or None, log_dict)
    """
    log = {
        "original_species": scientific_name,
        "mapped_species": "",
        "accession": "",
        "feature_location": "",
        "product": "",
        "sequence_length": "",
        "search_level": "GENUS_STANDIN",
        "status": "",
        "ncbi_url": "",
        "notes": "",
    }

    duplicates_seen: Set[str] = set()

    term = build_genus_search_term(scientific_name)
    if term is None:
        log["status"] = "NO_GENUS_TERM"
        log["notes"] = "Could not determine genus."
        return None, log

    print(f"[INFO] GENUS fallback term for '{scientific_name}': {term}")

    try:
        ids = entrez_search_nuccore(term, retmax=100)
    except Exception as e:
        log["status"] = "ERROR"
        log["notes"] = f"esearch genus failed: {e}"
        print(
            f"[ERROR] esearch genus failed for '{scientific_name}' "
            f"term '{term}': {e}"
        )
        return None, log

    if not ids:
        log["status"] = "NO_GENUS_HITS"
        log["notes"] = "No genus-level nuccore hits."
        return None, log

    for nid in ids:
        record = entrez_fetch_gb(nid)
        if record is None:
            continue

        extracted = extract_16s_feature(record)
        if extracted is None:
            continue

        seq_str, product, loc_str = extracted
        if seq_str in existing_seq_strings:
            existing_species = existing_seq_map.get(seq_str, set())
            duplicates_seen.update(existing_species)
            print(
                f"[INFO] GENUS candidate for '{scientific_name}' "
                f"is a duplicate of existing: {', '.join(existing_species)}; skipping."
            )
            continue

        mapped_species = record.annotations.get("organism", "").strip()
        if not mapped_species:
            mapped_species = scientific_name

        accession = record.id
        header = make_header(scientific_name, mapped_species)
        new_rec = SeqRecord(
            record.seq.__class__(seq_str),
            id=header,
            name=header,
            description=header,
        )

        existing_seq_strings.add(seq_str)

        log["mapped_species"] = mapped_species
        log["accession"] = accession
        log["feature_location"] = loc_str
        log["product"] = product
        log["sequence_length"] = str(len(seq_str))
        log["status"] = "GENUS_STANDIN_FOUND"
        log["ncbi_url"] = f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
        log["notes"] = ""
        return new_rec, log

    if duplicates_seen:
        log["status"] = "ONLY_DUPLICATE_16S_GENUS"
        log["notes"] = (
            "Genus-level search found only duplicate 16S sequences; "
            "existing species: " + ", ".join(sorted(duplicates_seen))
        )
    else:
        log["status"] = "NO_GENUS_16S"
        log["notes"] = (
            "No acceptable 16S found at genus level "
            "(no 16S features or no usable records)."
        )

    return None, log


# ------------------ MAIN ------------------ #

def main():
    existing_seq_strings, existing_seq_map = load_existing_sequences(COMBINED_FASTA)
    missing_species = load_missing_species(MISSING_TXT)

    rescued_records: List[SeqRecord] = []
    logs: List[dict] = []
    no_hits: List[str] = []

    for i, species in enumerate(missing_species, start=1):
        print(f"\n========== [{i}/{len(missing_species)}] {species} ==========")

        # 1. Try direct species
        rec, log = try_direct_species(species, existing_seq_strings, existing_seq_map)
        logs.append(log)

        if rec is not None:
            print(f"[SUCCESS] DIRECT 16S rescued for '{species}'.")
            rescued_records.append(rec)
            continue

        # 2. If no direct, try genus stand-in
        rec2, log2 = try_genus_standin(
            species, existing_seq_strings, existing_seq_map
        )
        logs.append(log2)

        if rec2 is not None:
            print(
                f"[SUCCESS] GENUS STAND-IN 16S rescued for '{species}' "
                f"using '{log2.get('mapped_species', '')}'."
            )
            rescued_records.append(rec2)
        else:
            status2 = log2.get("status", "")
            # If genus also only found duplicates, we still consider this a no-hit
            if status2.startswith("ONLY_DUPLICATE_16S"):
                print(
                    f"[FAIL] Only duplicate 16S sequences found for '{species}' "
                    f"at genus level."
                )
            else:
                print(f"[FAIL] No acceptable 16S found for '{species}'.")
            no_hits.append(species)

    # Write rescued FASTA
    if rescued_records:
        SeqIO.write(rescued_records, OUT_FASTA, "fasta")
        print(
            f"\n[INFO] Wrote {len(rescued_records)} rescued sequences to {OUT_FASTA}"
        )
    else:
        print("\n[INFO] No sequences were rescued; no FASTA written.")

    # Write full log CSV
    fieldnames = [
        "original_species",
        "mapped_species",
        "accession",
        "feature_location",
        "product",
        "sequence_length",
        "search_level",
        "status",
        "ncbi_url",
        "notes",
    ]
    with open(OUT_LOG, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in logs:
            for key in fieldnames:
                row.setdefault(key, "")
            writer.writerow(row)
    print(f"[INFO] Wrote log to {OUT_LOG}")

    # Write "no hits" list (species with no accepted new 16S)
    if no_hits:
        with open(OUT_NO_HITS, "w") as f:
            for name in no_hits:
                f.write(name + "\n")
        print(f"[INFO] Wrote {len(no_hits)} no-hit species to {OUT_NO_HITS}")
    else:
        print("[INFO] All species had some accepted 16S hit.")

    # Extract "duplicates-only" situations
    duplicates_only_rows = [
        row
        for row in logs
        if row.get("status", "") in (
            "ONLY_DUPLICATE_16S_DIRECT",
            "ONLY_DUPLICATE_16S_GENUS",
        )
    ]

    if duplicates_only_rows:
        with open(OUT_DUP_ONLY, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in duplicates_only_rows:
                for key in fieldnames:
                    row.setdefault(key, "")
                writer.writerow(row)
        print(
            f"[INFO] Wrote {len(duplicates_only_rows)} 'duplicates-only' "
            f"rows to {OUT_DUP_ONLY}"
        )
    else:
        print("[INFO] No 'duplicates-only' cases detected.")


if __name__ == "__main__":
    main()
