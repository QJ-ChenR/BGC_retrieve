"""Evidence-gated correction of legacy domain C-terminal coordinates.

This module plans changes without mutating records, then edits only two numeric
tokens in GenBank bytes. DNA, translations, other annotations and formatting are
preserved. It is deliberately not a general annotation repair engine.
"""

from __future__ import annotations

from hashlib import sha256
import json
import re

from Bio.SeqFeature import ExactPosition, SimpleLocation


RULE = "inclusive_protein_end_and_missing_terminal_codon"
GENOMIC_ONLY_RULE = "genomic_end_missing_terminal_codon"
IDENTIFIERS = ("locus_tag", "protein_id", "gene")
DOMAIN_TYPES = frozenset({"aSDomain", "PFAM_domain"})


def _digest(text):
    return sha256(text.encode("ascii")).hexdigest()


def _simple_exact(feature):
    loc = feature.location
    return (isinstance(loc, SimpleLocation) and len(loc.parts) == 1
            and type(loc.start) is ExactPosition and type(loc.end) is ExactPosition
            and loc.ref is None and loc.ref_db is None and loc.strand in (-1, 1))


def plan_repairs(record, protein_module):
    """Return (proven changes, checks for aSDomain/PFAM_domain), without mutation.

The current encoder's translation, location and parent-resolution checks are
reused. A repair additionally requires a unique geometric parent independently
of identifiers and matching annotated and DNA-derived peptide evidence.
"""
    parser = protein_module
    parents = []
    for index, feature in enumerate(record.features):
        if feature.type != "CDS":
            continue
        parents.append({
            **parser._location(feature, len(record.seq)),
            "feature_index": index, "feature": feature,
            "sequence": parser._translation(feature),
            "identifiers": {key: value for key in IDENTIFIERS
                            if (value := parser._one(feature, key)) is not None},
        })
    if not parents:
        raise ValueError("No CDS features found")
    changes, checks = [], []
    for index, feature in enumerate(record.features):
        if feature.type not in DOMAIN_TYPES:
            continue
        sequence = parser._translation(feature)
        location = parser._location(feature, len(record.seq))
        aliases = {value for key in IDENTIFIERS
                   if (value := parser._one(feature, key)) is not None}
        alias_matches = [p for p in parents if aliases.intersection(p["identifiers"].values())]
        geometric = [p for p in parents if parser._contained(location, p)]
        candidates = [p for p in (alias_matches or parents) if parser._contained(location, p)]
        if len(candidates) != 1:
            if feature.type == "PFAM_domain":
                checks.append({"feature_index": index, "feature_type": feature.type,
                               "status": "needs_review", "geometric_parent_count": len(geometric),
                               "reason": f"expected one parent, found {len(candidates)}; PFAM unchanged",
                               "domain_sha256": _digest(sequence), "domain_length": len(sequence)})
                continue
            raise ValueError(f"aSDomain {index}: expected one parent, found {len(candidates)}")
        parent = candidates[0]
        base = {
            "feature_index": index, "feature_type": feature.type,
            "parent_feature_index": parent["feature_index"],
            "domain_name": parser._one(feature, "aSDomain") or "",
            "parent_identifiers_json": json.dumps(parent["identifiers"], sort_keys=True, ensure_ascii=False),
            "parent_resolution": "identifiers_and_containment" if alias_matches else "containment",
            "geometric_parent_count": len(geometric), "domain_length": len(sequence),
            "domain_sha256": _digest(sequence),
            "parent_translation_sha256": _digest(parent["sequence"]),
            "strand": location["strand"],
            "old_genomic_start": location["genomic_start"],
            "old_genomic_end": location["genomic_end"],
        }
        try:
            start, end, source = parser._protein_interval(feature, sequence, parent["sequence"])
        except ValueError as original_error:
            start_text = parser._one(feature, "protein_start")
            end_text = parser._one(feature, "protein_end")
            if start_text is None or end_text is None:
                raise ValueError(f"aSDomain {index}: no explicit complete coordinates; {original_error}") from original_error
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as exc:
                raise ValueError(f"aSDomain {index}: coordinates are not integers") from exc
            resolved_end = end + 1
            if (not 0 <= start <= end < len(parent["sequence"])
                    or resolved_end - start != len(sequence)
                    or parent["sequence"][start:resolved_end] != sequence):
                raise ValueError(f"aSDomain {index}: end+1 lacks exact parent-sequence evidence") from original_error
            old_protein_end, new_protein_end = end, resolved_end
            source, rule = "annotation_inclusive_end_repaired", RULE
        else:
            resolved_end = end
            raw_length = len(feature.location)
            if raw_length == 3 * len(sequence):
                checks.append({**base, "status": "unchanged", "protein_start": start,
                               "old_protein_end": parser._one(feature, "protein_end"),
                               "new_protein_end": parser._one(feature, "protein_end"),
                               "resolved_protein_start": start, "resolved_protein_end": end,
                               "new_genomic_start": location["genomic_start"],
                               "new_genomic_end": location["genomic_end"],
                               "coordinate_source": source, "parent_slice_matches": True, "rule": ""})
                continue
            if source != "unique_translation_match":
                raise ValueError(f"aSDomain {index}: unexpected genomic discrepancy with valid explicit protein coordinates")
            old_protein_end = new_protein_end = None
            rule = GENOMIC_ONLY_RULE
        if len(geometric) != 1 or geometric[0]["feature_index"] != parent["feature_index"]:
            raise ValueError(f"aSDomain {index}: repair requires a unique geometric CDS parent")
        parent_feature = parent["feature"]
        if not _simple_exact(feature) or not _simple_exact(parent_feature):
            raise ValueError(f"aSDomain {index}: repair requires simple exact domain and CDS locations")
        strand = location["strand"]
        old_start, old_end = location["genomic_start"], location["genomic_end"]
        new_start, new_genomic_end = (old_start, old_end + 3) if strand == 1 else (old_start - 3, old_end)
        if not 0 <= new_start < new_genomic_end <= len(record.seq):
            raise ValueError(f"aSDomain {index}: proposed genomic extension leaves the record")
        proposed = SimpleLocation(new_start, new_genomic_end, strand=strand)
        proposed_feature = type(feature)(proposed, type="aSDomain")
        new_location = parser._location(proposed_feature, len(record.seq))
        if not parser._contained(new_location, parent):
            raise ValueError(f"aSDomain {index}: proposed genomic extension leaves parent CDS")
        try:
            codon_start = int(parser._one(parent_feature, "codon_start") or "1")
            table = int(parser._one(parent_feature, "transl_table") or "1")
        except ValueError as exc:
            raise ValueError(f"aSDomain {index}: invalid parent translation metadata") from exc
        if codon_start not in (1, 2, 3):
            raise ValueError(f"aSDomain {index}: invalid codon_start")
        offset = codon_start - 1
        if strand == 1:
            expected = (parent["genomic_start"] + offset + 3 * start,
                        parent["genomic_start"] + offset + 3 * resolved_end)
        else:
            expected = (parent["genomic_end"] - offset - 3 * resolved_end,
                        parent["genomic_end"] - offset - 3 * start)
        if (new_start, new_genomic_end) != expected:
            raise ValueError(f"aSDomain {index}: genomic extension disagrees with parent protein coordinates")
        if old_end - old_start != 3 * (len(sequence) - 1):
            raise ValueError(f"aSDomain {index}: original genomic length is not one codon short")
        original_peptide = str(feature.extract(record.seq).translate(table=table))
        repaired_peptide = str(proposed.extract(record.seq).translate(table=table))
        if original_peptide != sequence[:-1] or repaired_peptide != sequence:
            raise ValueError(f"aSDomain {index}: exact genomic translation evidence is absent")
        change = {
            **base, "protein_start": start, "old_protein_end": old_protein_end,
            "new_protein_end": new_protein_end,
            "resolved_protein_start": start, "resolved_protein_end": resolved_end,
            "coordinate_source": source,
            "new_genomic_start": new_start, "new_genomic_end": new_genomic_end,
            "codon_start": codon_start, "transl_table": table,
            "original_dna_translation_sha256": _digest(original_peptide),
            "repaired_dna_translation_sha256": _digest(repaired_peptide),
            "parent_slice_matches": True, "repaired_dna_translation_matches": True,
            "rule": rule,
        }
        changes.append(change)
        checks.append({**change, "status": "repaired"})
    return changes, checks


def apply_byte_repairs(raw_bytes, record, changes):
    """Change only proven numeric coordinate tokens, preserving all other bytes.

Record feature indices must exactly match the raw FEATURES section. Unsupported
multiline/simple-location formatting is rejected instead of reserializing GBK.
"""
    if not isinstance(raw_bytes, bytes):
        raise TypeError("raw_bytes must be bytes")
    sections = list(re.finditer(rb"(?m)^FEATURES[^\r\n]*\r?$", raw_bytes))
    origins = list(re.finditer(rb"(?m)^ORIGIN[^\r\n]*\r?$", raw_bytes))
    if len(sections) != 1 or len(origins) != 1 or sections[0].end() >= origins[0].start():
        raise ValueError("Expected one FEATURES section followed by one ORIGIN")
    section_start, section_end = sections[0].end(), origins[0].start()
    section = raw_bytes[section_start:section_end]
    headers = list(re.finditer(rb"(?m)^ {5}(\S+)[ \t]+([^\r\n]*)\r?$", section))
    if len(headers) != len(record.features):
        raise ValueError("Raw and parsed feature counts differ")
    if any(match[1].decode("ascii") != feature.type for match, feature in zip(headers, record.features)):
        raise ValueError("Raw and parsed feature types/order differ")
    seen, replacements = set(), []
    for change in changes:
        index = int(change["feature_index"])
        if index in seen or not 0 <= index < len(headers):
            raise ValueError("Duplicate or invalid repair feature index")
        seen.add(index)
        feature = record.features[index]
        if feature.type not in DOMAIN_TYPES or not _simple_exact(feature):
            raise ValueError("Only simple exact domain locations may be patched")
        if (int(feature.location.start), int(feature.location.end), feature.location.strand) != (
                change["old_genomic_start"], change["old_genomic_end"], change["strand"]):
            raise ValueError("Repair location does not match parsed source")
        if change.get("rule") == RULE:
            if feature.qualifiers.get("protein_end") != [str(change["old_protein_end"])]:
                raise ValueError("Repair protein_end does not match parsed source")
            if change["new_protein_end"] != change["old_protein_end"] + 1:
                raise ValueError("Unrecognized protein-end operation")
        elif change.get("rule") == GENOMIC_ONLY_RULE:
            if (feature.qualifiers.get("protein_start") or feature.qualifiers.get("protein_end")
                    or change["old_protein_end"] is not None or change["new_protein_end"] is not None):
                raise ValueError("Genomic-only operation must not invent protein qualifiers")
        else:
            raise ValueError("Unrecognized repair operation")
        header = headers[index]
        loc_text = header[2].rstrip(b" \t")
        positive = change["strand"] == 1
        pattern = rb"([0-9]+)\.\.([0-9]+)" if positive else rb"complement\(([0-9]+)\.\.([0-9]+)\)"
        loc_match = re.fullmatch(pattern, loc_text)
        if loc_match is None:
            raise ValueError("Unsupported aSDomain location text; refusing serialization")
        if (int(loc_match[1]) - 1, int(loc_match[2])) != (change["old_genomic_start"], change["old_genomic_end"]):
            raise ValueError("Raw location numeric tokens differ from parsed source")
        expected_new = ((change["old_genomic_start"], change["old_genomic_end"] + 3) if positive
                        else (change["old_genomic_start"] - 3, change["old_genomic_end"]))
        if (change["new_genomic_start"], change["new_genomic_end"]) != expected_new:
            raise ValueError("Unrecognized genomic extension operation")
        group = 2 if positive else 1
        start, end = loc_match.span(group)
        start += section_start + header.start(2)
        end += section_start + header.start(2)
        new_value = change["new_genomic_end"] if positive else change["new_genomic_start"] + 1
        replacements.append((start, end, str(new_value).encode("ascii")))
        change["source_location_line"] = raw_bytes.count(b"\n", 0, start) + 1
        if change["rule"] == GENOMIC_ONLY_RULE:
            change["source_protein_end_line"] = None
            continue
        block_start = section_start + header.end()
        block_end = section_start + (headers[index + 1].start() if index + 1 < len(headers) else len(section))
        block = raw_bytes[block_start:block_end]
        qualifiers = list(re.finditer(
            rb'(?m)^ {21}/protein_end=(?P<quote>"?)(?P<value>[0-9]+)(?P=quote)[ \t]*\r?$', block))
        if len(qualifiers) != 1 or int(qualifiers[0]["value"]) != change["old_protein_end"]:
            raise ValueError("Expected exactly one matching raw protein_end qualifier")
        value_start, value_end = qualifiers[0].span("value")
        value_start += block_start
        value_end += block_start
        replacements.append((value_start, value_end, str(change["new_protein_end"]).encode("ascii")))
        change["source_location_line"] = raw_bytes.count(b"\n", 0, start) + 1
        change["source_protein_end_line"] = raw_bytes.count(b"\n", 0, value_start) + 1
    replacements.sort()
    if any(left[1] > right[0] for left, right in zip(replacements, replacements[1:])):
        raise ValueError("Overlapping raw byte repairs")
    result = raw_bytes
    for start, end, value in reversed(replacements):
        result = result[:start] + value + result[end:]
    return result
