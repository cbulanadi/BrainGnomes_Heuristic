# -*- coding: utf-8 -*-
"""
CogMAP heuristic for HeuDiConv (session-safe, stricter matching)
"""

import re
import json
from pathlib import Path

# Auto-populate IntendedFor for fmaps using closest match and acquisition labels
POPULATE_INTENDED_FOR_OPTS = {
    "matching_parameters": ["CustomAcquisitionLabel"],
    "criterion": "Closest",
}


def create_key(template, outtype=("nii.gz",), annotation_classes=None):
    if not template:
        raise ValueError("Template must be a valid format string")
    return template, outtype, annotation_classes


def _task_run_match(desc, task, run_num):
    # Matches underscore/hyphen/space separated labels like:
    # ccf_run1, ccf-run-01, ccf part1, ccf_part_01, ahcp_fmri_dpx_run1, rise_part2
    task_present = re.search(rf"(^|[^a-z0-9]){task}($|[^a-z0-9])", desc)
    run_present = re.search(rf"(^|[^a-z0-9])(run|part)[^a-z0-9]*0?{run_num}($|[^0-9])", desc)
    return bool(task_present and run_present)


def _as_lower_text(value):
    if not value:
        return ""
    if isinstance(value, (tuple, list)):
        return " ".join(str(v) for v in value).lower()
    return str(value).lower()


def _is_derived_or_secondary(s):
    """
    Avoid assigning non-primary or derived images (which can trigger conversion
    failures like missing PixelSpacing and duplicate output collisions).
    """
    image_type = _as_lower_text(getattr(s, "image_type", ""))
    bad_tokens = ("derived", "secondary", "localizer", "scout", "moco")
    return any(tok in image_type for tok in bad_tokens)


def _is_plausible_image_series(s):
    """Basic sanity filter for real image volumes."""
    if _is_derived_or_secondary(s):
        return False

    dim1 = int(getattr(s, "dim1", 0) or 0)
    dim2 = int(getattr(s, "dim2", 0) or 0)
    return dim1 >= 32 and dim2 >= 32


def _is_probable_dwi(desc, s):
    """
    Keep only true diffusion series for *_dwi outputs to avoid
    VOLUME_COUNT_MISMATCH from assigning non-diffusion PA/AP scans.
    """
    full_text = " ".join(
        [
            desc,
            _as_lower_text(getattr(s, "protocol_name", "")),
            _as_lower_text(getattr(s, "sequence_name", "")),
            _as_lower_text(getattr(s, "image_type", "")),
        ]
    )

    excluded_terms = (
        "distortionmap",
        "fieldmap",
        "topup",
        "fmap",
        "sbref",
        "trace",
        "adc",
        "fa",
        "md",
    )
    if any(term in full_text for term in excluded_terms):
        return False

    dim4 = int(getattr(s, "dim4", 0) or 0)
    # True dMRI should have many volumes; avoid sending short AP/PA map scans to dwi.
    if dim4 < 10:
        return False

    return ("dmri" in full_text or "dwi" in full_text) and _is_plausible_image_series(s)


def _append_once(info, key, series_id):
    """
    Prevent duplicate writes to fixed output names (for example T1w, run-01,
    run-02), which can otherwise fail with destination-already-exists.
    """
    if not info[key]:
        info[key].append(series_id)


def _sanitize_bids_label(value):
    """Convert arbitrary text into a safe BIDS entity label fragment."""
    text = _as_lower_text(value)
    cleaned = re.sub(r"[^a-z0-9]+", "", text)
    return cleaned


def infotoids(seqinfos, outdir=None):
    """
    Populate subject/session/locator IDs from sequence metadata when available.

    This suppresses heudiconv's fallback warning about missing ``infotoids`` and
    ensures conversions remain stable across all subject IDs even when command
    line ``-s/-ss`` values are not explicitly provided.
    """
    if not seqinfos:
        return {"locator": "cogmap", "session": "ses-01"}

    first = seqinfos[0]

    # Subject: prefer explicit DICOM fields, then fallback to first numeric token.
    subject_raw = (
        getattr(first, "patient_id", None)
        or getattr(first, "patient_name", None)
        or getattr(first, "study_description", None)
        or ""
    )
    subject_label = _sanitize_bids_label(subject_raw)
    if not subject_label:
        description = _as_lower_text(getattr(first, "series_description", ""))
        match = re.search(r"\b(\d{3,})\b", description)
        if match:
            subject_label = match.group(1)

    # Session: keep existing ses-* tags if present, otherwise infer a stable default.
    session_raw = (
        getattr(first, "study_id", None)
        or getattr(first, "study_description", None)
        or ""
    )
    session_label = _sanitize_bids_label(session_raw)
    if session_label.startswith("ses") and len(session_label) > 3:
        session = f"ses-{session_label[3:]}"
    elif session_label:
        session = f"ses-{session_label}"
    else:
        session = "ses-01"

    # Locator groups related datasets under a stable project-level namespace.
    locator_raw = getattr(first, "study_description", None) or "cogmap"
    locator = _sanitize_bids_label(locator_raw) or "cogmap"

    ids = {"locator": locator, "session": session}
    if subject_label:
        ids["subject"] = subject_label
    return ids


def _is_probable_bold(desc, s):
    """
    Restrict task assignments to true fMRI time series.

    This avoids assigning SBRef/localizer/derived/map series to *_bold names,
    which can create duplicate outputs like *_bold1.nii.gz that are not
    BIDS-compliant.
    """
    full_text = " ".join(
        [
            desc,
            _as_lower_text(getattr(s, "protocol_name", "")),
            _as_lower_text(getattr(s, "sequence_name", "")),
            _as_lower_text(getattr(s, "image_type", "")),
        ]
    )

    excluded_terms = (
        "sbref",
        "distortionmap",
        "fieldmap",
        "fmap",
        "topup",
        "dwi",
        "dmri",
        "localizer",
        "scout",
    )
    if any(term in full_text for term in excluded_terms):
        return False

    dim4 = int(getattr(s, "dim4", 0) or 0)
    if dim4 < 20:
        return False

    return _is_plausible_image_series(s)


def tuneup_bids_json_files(json_files):
    """
    Post-process sidecars after conversion.

    Some source DICOMs include both RepetitionTime and AcquisitionDuration,
    which triggers the BIDS validator error
    REPETITION_TIME_AND_ACQUISITION_DURATION_MUTUALLY_EXCLUSIVE.
    Keep RepetitionTime and remove AcquisitionDuration whenever both are
    present in the same sidecar.

    HeuDiConv may pass only newly-created JSONs to this hook. When rerunning
    partial conversions (for example, a single session), older sidecars under
    the same session might still carry AcquisitionDuration. To keep reruns
    stable, also sweep sibling JSON sidecars in the same session tree.
    """

    def _remove_non_bids_duplicate(path):
        """Delete accidental fallback files like *_bold1(.json/.nii.gz)."""
        if not re.search(r"_bold\d+\.json$", path.name):
            return False

        stem = path.name[:-5]  # drop .json
        nii = path.with_name(f"{stem}.nii.gz")
        for extra in (path, nii):
            if extra.exists():
                extra.unlink()
        return True

    def _cleanup_sidecar(path):
        if not path.exists():
            return
        if _remove_non_bids_duplicate(path):
            return

        # Remove accidental non-BIDS duplicate functional files (e.g. *_bold1)
        # that can appear when conversion encounters pre-existing targets.
        if re.search(r"_bold\d+\.json$", path.name):
            for extra in (path, path.with_suffix(""), path.with_suffix("").with_suffix(".nii.gz")):
                if extra.exists():
                    extra.unlink()
            continue

        with path.open("r", encoding="utf-8") as fobj:
            sidecar = json.load(fobj)

        if "RepetitionTime" in sidecar and "AcquisitionDuration" in sidecar:
            sidecar.pop("AcquisitionDuration", None)
            with path.open("w", encoding="utf-8") as fobj:
                json.dump(sidecar, fobj, indent=2, sort_keys=True)
                fobj.write("\n")

    # De-duplicate while preserving order.
    direct_paths = []
    seen = set()
    for json_file in json_files:
        path = Path(json_file)
        if path in seen:
            continue
        seen.add(path)
        direct_paths.append(path)

    for path in direct_paths:
        _cleanup_sidecar(path)

    # Also sanitize all JSON sidecars in the same sub-*/ses-* tree so reruns
    # don't leave older files with mutually-exclusive metadata.
    session_roots = set()
    for path in direct_paths:
        parts = path.parts
        for idx, part in enumerate(parts):
            if part.startswith("ses-") and idx > 0 and parts[idx - 1].startswith("sub-"):
                session_roots.add(Path(*parts[: idx + 1]))
                break

    for root in sorted(session_roots):
        if not root.exists() or not root.is_dir():
            continue
        for sidecar in root.rglob("*.json"):
            _cleanup_sidecar(sidecar)
