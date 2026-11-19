"""Restore explicitly marked Geneious aSDomain features from sequence evidence.

Coordinates are never repaired here. A candidate needs a unique, same-strand
CDS, exact complete codons, and agreement between genomic translation and the
annotated parent peptide. Repeated annotations remain repeated annotations.
"""

from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from hashlib import sha256
from io import StringIO
import re

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition, SimpleLocation


MARKER = "Geneious type: aSDomain"
RULE = "geneious_asdomain_exact_codon_parent_slice"
IDENTIFIERS = ("locus_tag", "protein_id", "gene")
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZUO")


def _digest(sequence):
    return sha256(sequence.encode("ascii")).hexdigest()


def _one(feature, key, default=None):
    values = feature.qualifiers.get(key)
    if values is None:
        return default
    if not isinstance(values, (list, tuple)) or len(values) != 1 or not isinstance(values[0], str):
        raise ValueError(f"invalid or multiple {key} qualifiers")
    return values[0]


def _simple_exact(feature):
    location = feature.location
    return (isinstance(location, SimpleLocation) and len(location.parts) == 1
            and type(location.start) is ExactPosition and type(location.end) is ExactPosition
            and location.ref is None and location.ref_db is None
            and location.strand in (-1, 1))


def _contains(parent, child):
    """Conservative containment, including complex CDSs as ambiguous candidates."""
    if parent.location is None or child.location is None:
        return False
    return (parent.location.strand == child.location.strand
            and int(parent.location.start) <= int(child.location.start)
            and int(child.location.end) <= int(parent.location.end))


def _identifiers(feature):
    return {key: value for key in IDENTIFIERS if (value := _one(feature, key))}


def _candidate(feature):
    notes = feature.qualifiers.get("note", [])
    return (feature.type == "misc_feature" and isinstance(notes, (list, tuple))
            and any(note == MARKER for note in notes))


def plan_reconstruction(record):
    """Return JSON-serializable ``(changes, checks)`` without mutating the record.

    There is one check per exact Geneious marker. Status is ``restorable`` or
    ``needs_review``. Duplicate groups refer to candidate annotations sharing
    a resolved parent and genomic interval; they are reported, never collapsed.
    """
    parents = [(index, feature) for index, feature in enumerate(record.features)
               if feature.type == "CDS"]
    changes, checks = [], []
    groups = defaultdict(list)
    for index, feature in enumerate(record.features):
        if not _candidate(feature):
            continue
        check = {
            "feature_index": index, "feature_type": feature.type,
            "rule": RULE, "genomic_location": str(feature.location),
            "old_qualifiers": deepcopy(feature.qualifiers),
            "status": "needs_review", "reason": "",
            "duplicate_group": "", "duplicate_group_size": 1,
        }
        checks.append(check)
        try:
            if not _simple_exact(feature):
                raise ValueError("domain location is not simple, exact and stranded")
            start, end = int(feature.location.start), int(feature.location.end)
            strand = feature.location.strand
            check.update(genomic_start=start, genomic_end=end, strand=strand,
                         genomic_length=end - start)
            if not 0 <= start < end <= len(record.seq):
                raise ValueError("domain location lies outside the record or is empty")
            matches = [(p_index, parent) for p_index, parent in parents if _contains(parent, feature)]
            check["geometric_parent_count"] = len(matches)
            if len(matches) != 1:
                raise ValueError(f"expected one same-strand geometric CDS parent, found {len(matches)}")
            parent_index, parent = matches[0]
            check["parent_feature_index"] = parent_index
            groups[(parent_index, start, end, strand)].append(check)
            if not _simple_exact(parent):
                raise ValueError("parent location is not simple, exact and stranded")
            parent_start, parent_end = int(parent.location.start), int(parent.location.end)
            if not 0 <= parent_start < parent_end <= len(record.seq):
                raise ValueError("parent location lies outside the record")
            ids = _identifiers(parent)
            old_ids = _identifiers(feature)
            check.update(parent_identifiers=ids, old_identifiers=old_ids,
                         parent_genomic_start=parent_start, parent_genomic_end=parent_end)
            parent_translation = "".join((_one(parent, "translation") or "").split()).upper()
            if not parent_translation or set(parent_translation) - AMINO_ACIDS:
                raise ValueError("parent translation is missing or contains unsupported amino acids")
            check["parent_translation_sha256"] = _digest(parent_translation)
            codon_start = int(_one(parent, "codon_start", "1"))
            transl_table = int(_one(parent, "transl_table", "1"))
            if codon_start not in (1, 2, 3):
                raise ValueError("parent codon_start must be 1, 2 or 3")
            check.update(codon_start=codon_start, transl_table=transl_table)
            offset = ((start - parent_start) if strand == 1 else (parent_end - end)) - (codon_start - 1)
            check.update(coding_offset_nt=offset, coding_offset_remainder=offset % 3,
                         genomic_length_remainder=(end - start) % 3)
            if offset < 0 or offset % 3 or (end - start) % 3:
                raise ValueError("domain boundaries are not complete codons in the parent reading frame")
            protein_start = offset // 3
            protein_end = protein_start + (end - start) // 3
            check.update(protein_start=protein_start, protein_end=protein_end)
            if not 0 <= protein_start < protein_end <= len(parent_translation):
                raise ValueError("mapped protein interval lies outside the parent translation")
            peptide = str(feature.extract(record.seq).translate(table=transl_table)).upper()
            check.update(domain_length=len(peptide), translation_sha256=_digest(peptide),
                         dna_translation_sha256=_digest(peptide))
            if set(peptide) - AMINO_ACIDS:
                raise ValueError("domain genomic translation contains unsupported amino acids or stop codons")
            if peptide != parent_translation[protein_start:protein_end]:
                raise ValueError("domain genomic translation does not equal the parent peptide slice")
            check.update(parent_slice_matches=True, dna_translation_matches=True)
            expected = {"translation": [peptide], "protein_start": [str(protein_start)],
                        "protein_end": [str(protein_end)]}
            for key, values in expected.items():
                if key not in feature.qualifiers:
                    continue
                existing = _one(feature, key)
                normalized = "".join(existing.split()).upper() if key == "translation" else existing
                if normalized != values[0]:
                    raise ValueError(f"existing {key} disagrees with sequence-derived evidence")
            new_qualifiers = deepcopy(feature.qualifiers)
            new_qualifiers.update(expected)
            # H's exported domains retain locus_tag. Keep the original evidence
            # in checks and map to the actual parent's existing tag if necessary.
            # Other contradictory identifier types need review, not invention.
            locus_tag = ids.get("locus_tag")
            if old_ids.get("locus_tag") != locus_tag:
                if not locus_tag:
                    raise ValueError("cannot map locus_tag: the parent has no existing locus_tag")
                new_qualifiers["locus_tag"] = [locus_tag]
            for key in ("protein_id", "gene"):
                if key in old_ids and old_ids[key] != ids.get(key):
                    raise ValueError(f"conflicting domain {key} requires manual review")
            identifier_mapping = {
                "old_locus_tag": old_ids.get("locus_tag"),
                "resolved_locus_tag": locus_tag,
                "changed": old_ids.get("locus_tag") != locus_tag,
                "evidence": "unique_same_strand_containment_and_exact_dna_parent_slice",
            }
            check.update(status="restorable", new_feature_type="aSDomain",
                         new_qualifiers=new_qualifiers, identifier_mapping=identifier_mapping,
                         parent_resolution="unique_same_strand_geometry_then_sequence",
                         translation=peptide)
            changes.append(check)
        except (ValueError, TypeError, KeyError) as error:
            check["reason"] = str(error)
    for key, members in groups.items():
        if len(members) > 1:
            group = "parent_{}_{}_{}_{}".format(*key)
            for check in members:
                check.update(duplicate_group=group, duplicate_group_size=len(members))
    return deepcopy(changes), checks


def _feature_spans(raw_bytes):
    sections = list(re.finditer(rb"(?m)^FEATURES[^\r\n]*\r?$", raw_bytes))
    origins = list(re.finditer(rb"(?m)^ORIGIN[^\r\n]*\r?$", raw_bytes))
    if len(sections) != 1 or len(origins) != 1 or sections[0].end() >= origins[0].start():
        raise ValueError("Expected one FEATURES section followed by one ORIGIN")
    start, end = sections[0].end(), origins[0].start()
    headers = list(re.finditer(rb"(?m)^ {5}(\S+)[ \t]+([^\r\n]*)\r?$", raw_bytes[start:end]))
    return [(start + header.start(), start + headers[i + 1].start() if i + 1 < len(headers) else end,
             header[1].decode("ascii")) for i, header in enumerate(headers)]


def _qualifier_bytes(key, value, newline):
    indent = b" " * 21
    if key in ("protein_start", "protein_end"):
        return indent + f"/{key}={value}".encode("ascii") + newline
    if key == "translation":
        prefix = indent + b'/translation="'
        chunks = [prefix + value[:44].encode("ascii")]
        rest = value[44:]
        chunks.extend(indent + rest[i:i + 58].encode("ascii") for i in range(0, len(rest), 58))
        chunks[-1] += b'"'
        return newline.join(chunks) + newline
    value = value.replace('"', '""')
    if "\n" in value or "\r" in value:
        raise ValueError("New locus_tag contains a newline")
    return indent + f'/{key}="{value}"'.encode("ascii") + newline


def _patched_feature(block, change):
    header = re.match(rb" {5}(\S+)([ \t]+)([^\r\n]*)(\r?\n)", block)
    if header is None or header[1] != b"misc_feature":
        raise ValueError("Unsupported source feature header")
    newline = header[4]
    header_bytes = b" " * 5 + b"aSDomain".ljust(len(header[1]) + len(header[2])) + header[3] + newline
    body = block[header.end():]
    starts = list(re.finditer(rb"(?m)^ {21}/([A-Za-z0-9_]+)(?:=|[ \t]*\r?$)", body))
    old = change["old_qualifiers"]
    new = change["new_qualifiers"]
    fields = [key for key in ("locus_tag", "translation", "protein_start", "protein_end")
              if old.get(key) != new.get(key)]
    replacements = []
    found = set()
    for i, match in enumerate(starts):
        key = match[1].decode("ascii")
        if key not in fields:
            continue
        if key in found or len(new[key]) != 1:
            raise ValueError("Repeated target qualifier cannot be safely patched")
        found.add(key)
        end = starts[i + 1].start() if i + 1 < len(starts) else len(body)
        replacements.append((match.start(), end, _qualifier_bytes(key, new[key][0], newline)))
    for start, end, replacement in reversed(replacements):
        body = body[:start] + replacement + body[end:]
    for key in fields:
        if key not in found:
            if key in old:
                raise ValueError(f"Parsed {key} has no matching raw qualifier block")
            body += _qualifier_bytes(key, new[key][0], newline)
    return header_bytes + body


def _same_record(before, after, changes):
    for key in ("id", "name", "description", "annotations", "dbxrefs", "letter_annotations"):
        if getattr(before, key) != getattr(after, key):
            raise ValueError(f"Record {key} changed during reconstruction")
    if str(before.seq) != str(after.seq) or len(before.features) != len(after.features):
        raise ValueError("Sequence or feature count changed during reconstruction")
    expected = {change["feature_index"]: change for change in changes}
    for index, (old, new) in enumerate(zip(before.features, after.features)):
        if old.location != new.location or old.id != new.id:
            raise ValueError(f"Feature {index} location or identity changed")
        change = expected.get(index)
        if change is None:
            expected_type, expected_qualifiers = old.type, old.qualifiers
        else:
            expected_type, expected_qualifiers = change["new_feature_type"], change["new_qualifiers"]
        if new.type != expected_type or new.qualifiers != expected_qualifiers:
            raise ValueError(f"Feature {index} differs from the planned reconstruction")


def apply_reconstruction_bytes(raw_bytes, record, changes):
    """Apply an evidence-checked plan, preserving all other GenBank bytes.

    The supplied record must actually describe these bytes. Plans are recomputed
    before patching, preventing stale or manually altered peptide/ID changes.
    The result is reparsed and checked against every feature and record field.
    """
    if not isinstance(raw_bytes, bytes):
        raise TypeError("raw_bytes must be bytes")
    parsed = SeqIO.read(StringIO(raw_bytes.decode("utf-8")), "genbank")
    _same_record(record, parsed, [])
    approved = {change["feature_index"]: change for change in plan_reconstruction(record)[0]}
    seen = set()
    for change in changes:
        index = change["feature_index"]
        if index in seen or approved.get(index) != change:
            raise ValueError("Duplicate, stale or unproven reconstruction change")
        seen.add(index)
    if not changes:
        return raw_bytes
    spans = _feature_spans(raw_bytes)
    if len(spans) != len(record.features) or any(
            feature.type != span[2] for feature, span in zip(record.features, spans)):
        raise ValueError("Raw and parsed feature counts/types/order differ")
    output = raw_bytes
    for change in sorted(changes, key=lambda item: item["feature_index"], reverse=True):
        start, end, _ = spans[change["feature_index"]]
        output = output[:start] + _patched_feature(output[start:end], change) + output[end:]
    restored = SeqIO.read(StringIO(output.decode("utf-8")), "genbank")
    _same_record(record, restored, changes)
    return output
