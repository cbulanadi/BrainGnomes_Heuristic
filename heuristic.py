"""Heuristic for BrainGnomes bids_conversion."""


def create_key(template, outtype=("nii.gz",), annotation_classes=None):
    if not template:
        raise ValueError("Template must be a valid format string")
    return template, outtype, annotation_classes


def infotodict(seqinfo):
    t1w = create_key("sub-{subject}/{session}/anat/sub-{subject}_{session}_T1w")

    task_ccf_run1 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-ccf_run-01_bold"
    )
    task_ccf_run1_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-ccf_run-01_sbref"
    )
    task_ccf_run2 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-ccf_run-02_bold"
    )
    task_ccf_run2_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-ccf_run-02_sbref"
    )
    fmap_ccf_pa = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-ccf_dir-PA_epi"
    )
    fmap_ccf_ap = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-ccf_dir-AP_epi"
    )

    task_dpx_run1 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-dpx_run-01_bold"
    )
    task_dpx_run1_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-dpx_run-01_sbref"
    )
    task_dpx_run2 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-dpx_run-02_bold"
    )
    task_dpx_run2_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-dpx_run-02_sbref"
    )
    fmap_dpx_pa = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dpx_dir-PA_epi"
    )
    fmap_dpx_ap = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dpx_dir-AP_epi"
    )

    task_rise_run1 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rise_run-01_bold"
    )
    task_rise_run1_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rise_run-01_sbref"
    )
    task_rise_run2 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rise_run-02_bold"
    )
    task_rise_run2_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rise_run-02_sbref"
    )
    fmap_rise_pa = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rise_dir-PA_epi"
    )
    fmap_rise_ap = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rise_dir-AP_epi"
    )

    task_emo_run1 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-emo_run-01_bold"
    )
    task_emo_run1_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-emo_run-01_sbref"
    )
    task_emo_run2 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-emo_run-02_bold"
    )
    task_emo_run2_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-emo_run-02_sbref"
    )
    fmap_emo_pa = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-emo_dir-PA_epi"
    )
    fmap_emo_ap = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-emo_dir-AP_epi"
    )

    task_rest_run1 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-01_bold"
    )
    task_rest_run1_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-01_sbref"
    )
    task_rest_run2 = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-02_bold"
    )
    task_rest_run2_sbref = create_key(
        "sub-{subject}/{session}/func/sub-{subject}_{session}_task-rest_run-02_sbref"
    )
    fmap_rest_pa = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rest_dir-PA_epi"
    )
    fmap_rest_ap = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-rest_dir-AP_epi"
    )

    dwi_ap = create_key("sub-{subject}/{session}/dwi/sub-{subject}_{session}_dir-AP_dwi")
    dwi_pa = create_key("sub-{subject}/{session}/dwi/sub-{subject}_{session}_dir-PA_dwi")
    dwi_sbref_ap = create_key(
        "sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-sbref_dir-AP_dwi"
    )
    dwi_sbref_pa = create_key(
        "sub-{subject}/{session}/dwi/sub-{subject}_{session}_acq-sbref_dir-PA_dwi"
    )
    fmap_dwi_ap = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dwi_dir-AP_epi"
    )
    fmap_dwi_pa = create_key(
        "sub-{subject}/{session}/fmap/sub-{subject}_{session}_acq-dwi_dir-PA_epi"
    )

    info = {
        t1w: [],
        task_ccf_run1: [],
        task_ccf_run1_sbref: [],
        task_ccf_run2: [],
        task_ccf_run2_sbref: [],
        fmap_ccf_pa: [],
        fmap_ccf_ap: [],
        task_dpx_run1: [],
        task_dpx_run1_sbref: [],
        task_dpx_run2: [],
        task_dpx_run2_sbref: [],
        fmap_dpx_pa: [],
        fmap_dpx_ap: [],
        task_rise_run1: [],
        task_rise_run1_sbref: [],
        task_rise_run2: [],
        task_rise_run2_sbref: [],
        fmap_rise_pa: [],
        fmap_rise_ap: [],
        task_emo_run1: [],
        task_emo_run1_sbref: [],
        task_emo_run2: [],
        task_emo_run2_sbref: [],
        fmap_emo_pa: [],
        fmap_emo_ap: [],
        task_rest_run1: [],
        task_rest_run1_sbref: [],
        task_rest_run2: [],
        task_rest_run2_sbref: [],
        fmap_rest_pa: [],
        fmap_rest_ap: [],
        dwi_ap: [],
        dwi_pa: [],
        dwi_sbref_ap: [],
        dwi_sbref_pa: [],
        fmap_dwi_ap: [],
        fmap_dwi_pa: [],
    }

    for s in seqinfo:
        desc = (s.series_description or "").lower()

        if "adni3_t1_mprag_sag_p2_iso" in desc:
            info[t1w].append(s.series_id)

        elif "fmri_ccf_run1" in desc and "sbref" in desc:
            info[task_ccf_run1_sbref].append(s.series_id)
        elif "fmri_ccf_run1" in desc:
            info[task_ccf_run1].append(s.series_id)
        elif "fmri_ccf_run2" in desc and "sbref" in desc:
            info[task_ccf_run2_sbref].append(s.series_id)
        elif "fmri_ccf_run2" in desc:
            info[task_ccf_run2].append(s.series_id)
        elif "fmri_distortionmap_pa_ccf" in desc:
            info[fmap_ccf_pa].append(s.series_id)
        elif "fmri_distortionmap_ap_ccf" in desc:
            info[fmap_ccf_ap].append(s.series_id)

        elif "fmri_dpx_run1" in desc and "sbref" in desc:
            info[task_dpx_run1_sbref].append(s.series_id)
        elif "fmri_dpx_run1" in desc:
            info[task_dpx_run1].append(s.series_id)
        elif "fmri_dpx_run2" in desc and "sbref" in desc:
            info[task_dpx_run2_sbref].append(s.series_id)
        elif "fmri_dpx_run2" in desc:
            info[task_dpx_run2].append(s.series_id)
        elif "fmri_distortionmap_pa_dpx" in desc:
            info[fmap_dpx_pa].append(s.series_id)
        elif "fmri_distortionmap_ap_dpx" in desc:
            info[fmap_dpx_ap].append(s.series_id)

        elif "fmri_rise_part1" in desc and "sbref" in desc:
            info[task_rise_run1_sbref].append(s.series_id)
        elif "fmri_rise_part1" in desc:
            info[task_rise_run1].append(s.series_id)
        elif "fmri_rise_part2" in desc and "sbref" in desc:
            info[task_rise_run2_sbref].append(s.series_id)
        elif "fmri_rise_part2" in desc:
            info[task_rise_run2].append(s.series_id)
        elif "fmri_distortionmap_pa_rise" in desc:
            info[fmap_rise_pa].append(s.series_id)
        elif "fmri_distortionmap_ap_rise" in desc:
            info[fmap_rise_ap].append(s.series_id)

        elif "fmri_emo_run1" in desc and "sbref" in desc:
            info[task_emo_run1_sbref].append(s.series_id)
        elif "fmri_emo_run1" in desc:
            info[task_emo_run1].append(s.series_id)
        elif "fmri_emo_run2" in desc and "sbref" in desc:
            info[task_emo_run2_sbref].append(s.series_id)
        elif "fmri_emo_run2" in desc:
            info[task_emo_run2].append(s.series_id)
        elif "fmri_distortionmap_pa_emo" in desc:
            info[fmap_emo_pa].append(s.series_id)
        elif "fmri_distortionmap_ap_emo" in desc:
            info[fmap_emo_ap].append(s.series_id)

        elif "fmri_rest_run1" in desc and "sbref" in desc:
            info[task_rest_run1_sbref].append(s.series_id)
        elif "fmri_rest_run1" in desc:
            info[task_rest_run1].append(s.series_id)
        elif "fmri_rest_run2" in desc and "sbref" in desc:
            info[task_rest_run2_sbref].append(s.series_id)
        elif "fmri_rest_run2" in desc:
            info[task_rest_run2].append(s.series_id)
        elif "fmri_distortionmap_pa_rest" in desc:
            info[fmap_rest_pa].append(s.series_id)
        elif "fmri_distortionmap_ap_rest" in desc:
            info[fmap_rest_ap].append(s.series_id)

        elif "dmri_distortionmap_pa" in desc:
            info[fmap_dwi_pa].append(s.series_id)
        elif "dmri_distortionmap_ap" in desc:
            info[fmap_dwi_ap].append(s.series_id)
        elif "dmri_pa_sbref" in desc:
            info[dwi_sbref_pa].append(s.series_id)
        elif "dmri_ap_sbref" in desc:
            info[dwi_sbref_ap].append(s.series_id)
        elif desc == "dmri_pa" or desc.endswith("_dmri_pa"):
            info[dwi_pa].append(s.series_id)
        elif desc == "dmri_ap" or desc.endswith("_dmri_ap"):
            info[dwi_ap].append(s.series_id)

    return info
