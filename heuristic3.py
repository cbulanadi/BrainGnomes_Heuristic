# -*- coding: utf-8 -*-
"""
CogMAP heuristic for HeuDiConv (session-safe, stricter matching)
"""

import re

# Auto-populate IntendedFor for fmaps using acquisition labels.
#
# NOTE:
# Some sessions can yield sidecars with AcquisitionTime="n/a" (for example when
# dcmstack embedding fails on missing PixelSpacing in source DICOM metadata).
# HeuDiConv's "Closest" criterion requires parsing AcquisitionTime and raises:
#   ValueError: Unable to parse datetime string: n/a
# Using "First" avoids datetime parsing while still assigning IntendedFor
# within compatible acquisition groups.
POPULATE_INTENDED_FOR_OPTS = {
    "matching_parameters": ["CustomAcquisitionLabel"],
    "criterion": "First",
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


def _full_text(desc, s):
    return " ".join(
        [
            desc,
            _as_lower_text(getattr(s, "protocol_name", "")),
            _as_lower_text(getattr(s, "sequence_name", "")),
            _as_lower_text(getattr(s, "image_type", "")),
        ]
    )


def _is_probable_dwi(desc, s):
    """
    Keep only true diffusion series for *_dwi outputs to avoid
    VOLUME_COUNT_MISMATCH from assigning non-diffusion PA/AP scans.
    """
    full_text = _full_text(desc, s)

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

    return "dmri" in full_text or "dwi" in full_text


def _is_probable_fmap(full_text):
    fmap_terms = ("distortionmap", "fieldmap", "topup", "fmap")
    return any(term in full_text for term in fmap_terms)


def infotoids(seqinfos, outdir):
    """
    Provide stable subject/session fields to avoid ambiguity in BIDS path building.
    """
    return {
        "locator": None,
        "session": None,
        "subject": None,
    }


def infotodict(seqinfo):
    # IMPORTANT:
    # Use {session} directly (do NOT prepend ses-), because in this project
    # session is already coming through as ses-01 / ses-02.
    t1w = create_key("sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w")

    task_ccf_1 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-ccf_run-01_bold")
    task_ccf_2 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-ccf_run-02_bold")
    task_rise_1 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-rise_run-01_bold")
    task_rise_2 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-rise_run-02_bold")
    task_dpx_1 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-dpx_run-01_bold")
    task_dpx_2 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-dpx_run-02_bold")
    task_emo_1 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-emo_run-01_bold")
    task_emo_2 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-emo_run-02_bold")
    task_rest_1 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-01_bold")
    task_rest_2 = create_key("sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-02_bold")

    # Include run-{item:02d} to prevent multiple AP/PA diffusion series from colliding.
    dwi_ap = create_key("sub-{subject}/{session}/dwi/sub-{subject}_{session}_dir-AP_run-{item:02d}_dwi")
    dwi_pa = create_key("sub-{subject}/{session}/dwi/sub-{subject}_{session}_dir-PA_run-{item:02d}_dwi")

    # Include run-{item:02d} on all fmap outputs to prevent collisions when
    # multiple AP/PA distortion maps share the same acquisition label.
    fmap_ap = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-AP_run-{item:02d}_epi")
    fmap_pa = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-PA_run-{item:02d}_epi")
    fmap_ap_ccf = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-ccf_dir-AP_run-{item:02d}_epi")
    fmap_pa_ccf = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-ccf_dir-PA_run-{item:02d}_epi")
    fmap_ap_rise = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rise_dir-AP_run-{item:02d}_epi")
    fmap_pa_rise = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rise_dir-PA_run-{item:02d}_epi")
    fmap_ap_dpx = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dpx_dir-AP_run-{item:02d}_epi")
    fmap_pa_dpx = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dpx_dir-PA_run-{item:02d}_epi")
    fmap_ap_emo = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-emo_dir-AP_run-{item:02d}_epi")
    fmap_pa_emo = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-emo_dir-PA_run-{item:02d}_epi")

    info = {
        t1w: [],
        task_ccf_1: [], task_ccf_2: [],
        task_rise_1: [], task_rise_2: [],
        task_dpx_1: [], task_dpx_2: [],
        task_emo_1: [], task_emo_2: [],
        task_rest_1: [], task_rest_2: [],
        dwi_ap: [], dwi_pa: [],
        fmap_ap: [], fmap_pa: [],
        fmap_ap_ccf: [], fmap_pa_ccf: [],
        fmap_ap_rise: [], fmap_pa_rise: [],
        fmap_ap_dpx: [], fmap_pa_dpx: [],
        fmap_ap_emo: [], fmap_pa_emo: [],
    }

    for s in seqinfo:
        desc = (s.series_description or "").lower()
        full_text = _full_text(desc, s)
        is_fmap = _is_probable_fmap(full_text)

        # Ignore localizers/scouts
        if "localizer" in full_text or "scout" in full_text:
            continue

        # Fieldmaps first (to avoid accidental task matches)
        if is_fmap and "ap" in full_text and "ccf" in full_text:
            info[fmap_ap_ccf].append(s.series_id)
        elif is_fmap and "pa" in full_text and "ccf" in full_text:
            info[fmap_pa_ccf].append(s.series_id)
        elif is_fmap and "ap" in full_text and "rise" in full_text:
            info[fmap_ap_rise].append(s.series_id)
        elif is_fmap and "pa" in full_text and "rise" in full_text:
            info[fmap_pa_rise].append(s.series_id)
        elif is_fmap and "ap" in full_text and "dpx" in full_text:
            info[fmap_ap_dpx].append(s.series_id)
        elif is_fmap and "pa" in full_text and "dpx" in full_text:
            info[fmap_pa_dpx].append(s.series_id)
        elif is_fmap and "ap" in full_text and "emo" in full_text:
            info[fmap_ap_emo].append(s.series_id)
        elif is_fmap and "pa" in full_text and "emo" in full_text:
            info[fmap_pa_emo].append(s.series_id)
        elif is_fmap and "ap" in full_text:
            info[fmap_ap].append(s.series_id)
        elif is_fmap and "pa" in full_text:
            info[fmap_pa].append(s.series_id)

        # Diffusion
        elif ("dmri_ap" in full_text or "dwi_ap" in full_text) and _is_probable_dwi(desc, s):
            info[dwi_ap].append(s.series_id)
        elif ("dmri_pa" in full_text or "dwi_pa" in full_text) and _is_probable_dwi(desc, s):
            info[dwi_pa].append(s.series_id)

        # Anatomical
        elif re.search(r"(^|[^a-z0-9])(t1|t1w|mpr|mprag|mprage)($|[^a-z0-9])", desc):
            info[t1w].append(s.series_id)

        # Functional tasks
        elif _task_run_match(desc, "ccf", 1):
            info[task_ccf_1].append(s.series_id)
        elif _task_run_match(desc, "ccf", 2):
            info[task_ccf_2].append(s.series_id)
        elif _task_run_match(desc, "rise", 1):
            info[task_rise_1].append(s.series_id)
        elif _task_run_match(desc, "rise", 2):
            info[task_rise_2].append(s.series_id)
        elif _task_run_match(desc, "dpx", 1):
            info[task_dpx_1].append(s.series_id)
        elif _task_run_match(desc, "dpx", 2):
            info[task_dpx_2].append(s.series_id)
        elif _task_run_match(desc, "emo", 1):
            info[task_emo_1].append(s.series_id)
        elif _task_run_match(desc, "emo", 2):
            info[task_emo_2].append(s.series_id)
        elif _task_run_match(desc, "rest", 1):
            info[task_rest_1].append(s.series_id)
        elif _task_run_match(desc, "rest", 2):
            info[task_rest_2].append(s.series_id)

    return info
