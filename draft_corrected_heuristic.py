# -*- coding: utf-8 -*-
"""
CogMAP heuristic for HeuDiConv (session-safe, stricter matching)
"""

import re

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
    # Matches examples like:
    # ccf_run1, ccf-run-01, ccf part1, ccf_part_01, etc.
    patterns = [
        rf"\b{task}[_\-\s]*(run|part)?[_\-\s]*0?{run_num}\b",
        rf"\b{task}\b.*\b(run|part)\b[_\-\s]*0?{run_num}\b",
    ]
    return any(re.search(p, desc) for p in patterns)


def _as_lower_text(value):
    if not value:
        return ""
    if isinstance(value, (tuple, list)):
        return " ".join(str(v) for v in value).lower()
    return str(value).lower()


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

    return "dmri" in full_text or "dwi" in full_text


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

    fmap_ap = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-AP_epi")
    fmap_pa = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_dir-PA_epi")
    fmap_ap_ccf = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-ccf_dir-AP_epi")
    fmap_pa_ccf = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-ccf_dir-PA_epi")
    fmap_ap_rise = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rise_dir-AP_epi")
    fmap_pa_rise = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rise_dir-PA_epi")
    fmap_ap_dpx = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dpx_dir-AP_epi")
    fmap_pa_dpx = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dpx_dir-PA_epi")
    fmap_ap_emo = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-emo_dir-AP_epi")
    fmap_pa_emo = create_key("sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-emo_dir-PA_epi")

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

        # Ignore localizers/scouts
        if "localizer" in desc or "scout" in desc:
            continue

        # Fieldmaps first (to avoid accidental task matches)
        if "distortionmap" in desc and "ap" in desc and "ccf" in desc:
            info[fmap_ap_ccf].append(s.series_id)
        elif "distortionmap" in desc and "pa" in desc and "ccf" in desc:
            info[fmap_pa_ccf].append(s.series_id)
        elif "distortionmap" in desc and "ap" in desc and "rise" in desc:
            info[fmap_ap_rise].append(s.series_id)
        elif "distortionmap" in desc and "pa" in desc and "rise" in desc:
            info[fmap_pa_rise].append(s.series_id)
        elif "distortionmap" in desc and "ap" in desc and "dpx" in desc:
            info[fmap_ap_dpx].append(s.series_id)
        elif "distortionmap" in desc and "pa" in desc and "dpx" in desc:
            info[fmap_pa_dpx].append(s.series_id)
        elif "distortionmap" in desc and "ap" in desc and "emo" in desc:
            info[fmap_ap_emo].append(s.series_id)
        elif "distortionmap" in desc and "pa" in desc and "emo" in desc:
            info[fmap_pa_emo].append(s.series_id)
        elif "distortionmap" in desc and "ap" in desc:
            info[fmap_ap].append(s.series_id)
        elif "distortionmap" in desc and "pa" in desc:
            info[fmap_pa].append(s.series_id)

        # Diffusion
        elif ("dmri_ap" in desc or "dwi_ap" in desc) and _is_probable_dwi(desc, s):
            info[dwi_ap].append(s.series_id)
        elif ("dmri_pa" in desc or "dwi_pa" in desc) and _is_probable_dwi(desc, s):
            info[dwi_pa].append(s.series_id)

        # Anatomical
        elif re.search(r"\b(t1|t1w|mpr|mprage)\b", desc):
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
