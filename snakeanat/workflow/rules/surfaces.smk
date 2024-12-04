def get_parallel_opts_fastsurfer_seg(wildcards, threads, resources):
    if resources.gpu >= 1:
        return ""
    else: 
        return f"--threads {threads}"

def get_parallel_opts_fastsurfer_surf(wildcards, threads):
    return f"--threads {int(threads/2)} --parallel"


rule fastsurfer_seg:
    input: 
        t1=rules.get_template_t1.output,
    output:
        directory(
            sourcedata / "fastsurfer" / Path(bids(**inputs.subj_wildcards)).name
        )
    benchmark: out/f"code/benchmark/fastsurfer/{uid}.tsv"
    log: out/f"code/log/fastsurfer/{uid}.log"
    resources:
        gpu=1 if config.get('use_gpu',False) else 0,
        runtime=15,
        mem_mb=10000,
    container: config["containers"]["fastsurfer"],
    params:
        fs_license=os.environ.get('FS_LICENSE',config["fs_license"]),
        sid=lambda wcards, output: Path(output[0]).name,
        sd=lambda wcards, output: Path(output[0]).parent,
        parallel_opts=get_parallel_opts_fastsurfer_seg,
    threads: 
        1 if config.get('use_gpu',False) else 32
    group: 'fastsurfer_seg'
    shell:
        """
        /fastsurfer/run_fastsurfer.sh --fs_license {params.fs_license} \\
            --t1 {input.t1} --sid {params.sid} --sd {params.sd} \\
            --seg_only {params.parallel_opts} &> {log}
        """

rule fastsurfer_surf:
    input:
        t1=rules.get_template_t1.output,
        seg=rules.fastsurfer_seg.output,
    output:
        fs_dir=directory(
            sourcedata / "fastsurfer_surf" / Path(bids(**inputs.subj_wildcards)).name
        )
    benchmark: out/f"code/benchmark/fastsurfer_surf/{uid}.tsv"
    log: out/f"code/log/fastsurfer_surf/{uid}.log"
    container: config["containers"]["fastsurfer"]
    resources:
        runtime=150,
        mem_mb=10000,
    params:
        fs_license=os.environ.get('FS_LICENSE',config["fs_license"]),
        sid=lambda wcards, output: Path(output[0]).name,
        sd=lambda wcards, output: Path(output[0]).parent,
        parallel_opts=get_parallel_opts_fastsurfer_surf,
    threads: 16
    group: 'fastsurfer_surf'
    shell:
        """
        src="{input.seg}"
        length=$(printf '%s' "$src" | wc -m)
        for f in $(find "$src" -type f -o -type l); do
            rel="${{f:length}}"
            dest="{output}/$rel"
            dir=$(dirname "$dest")
            mkdir -p "$dir"
            ln -s "$(realpath --relative-to "$dir" "$src/$rel")" "$dest"
        done
        /fastsurfer/run_fastsurfer.sh --fs_license {params.fs_license} \\
            --t1 {input.t1} --sid {params.sid} --sd {params.sd} \\
            --surf_only {params.parallel_opts} &> {log}
        """

rule ciftify:
    input:
        fs_dir=rules.fastsurfer_surf.output.fs_dir,
    output:
        directory(
            sourcedata / "ciftify" / Path(bids(**inputs.subj_wildcards)).name
        )
    benchmark: out/f"code/benchmark/ciftify/{uid}.tsv"
    log: out/f"code/log/ciftify/{uid}.log"
    resources:
        runtime=240,
        mem_mb=10000,
    threads: 4
    params:
        sid=lambda wcards, output: Path(output[0]).name,
        sd=lambda wcards, output: Path(output[0]).parent,
        fs_dir=lambda wcards, input: Path(input[0]).parent,
        fs_license=os.environ.get('FS_LICENSE',config["fs_license"]),
    group: 'ciftify'
    container: config["containers"]["ciftify"]
    shell:
        """
        ciftify_recon_all {params.sid} \\
            --ciftify-work-dir {params.sd} --fs-subjects-dir {params.fs_dir}  \\
            --fs-license {params.fs_license} --n_cpus {threads} --resample-to-T1w32k \\
            --debug &> {log}
            
        """


rule bidsify:
    input: 
        data=rules.ciftify.output,
        config=Path(workflow.basedir) / '..' / 'resources' / 'from-ciftify_to-bids.yaml'
    output:
        flag=touch(
            sourcedata / "bidsify" / Path(bids(**inputs.subj_wildcards)).name
        ),
        invwarp=bids(
            out,
            datatype="xfm",
            mode="image",
            from_="MNI152NLin6Asym",
            to="T1w",
            suffix="xfm.nii.gz",
            **inputs.subj_wildcards,
        ),
        warp=bids(
            out,
            datatype="xfm",
            mode="image",
            to="MNI152NLin6Asym",
            from_="T1w",
            suffix="xfm.nii.gz",
            **inputs.subj_wildcards,
        ),
        affine=bids(
            out,
            datatype="xfm",
            mode="image",
            from_="T1w",
            to="MNI152NLin6Asym",
            suffix="xfm.mat",
            **inputs.subj_wildcards,
        ),
        mni_ref=bids(
            out,
            datatype="anat",
            space="MNI152NLin6ASym",
            desc="preproc",
            suffix="T1w.nii.gz",
            **inputs.subj_wildcards,
        ),
        t1w_ref=bids(
            out,
            datatype="anat",
            space="T1w",
            desc="preproc",
            suffix="T1w.nii.gz",
            **inputs.subj_wildcards,
        ),
    params:
        entities=lambda wcards: " ".join([
            "--bids " + "=".join(pair)
            for pair in inputs.subj_wildcards.items()
        ]).format(**wcards)
    container: config['containers']['snakeanat']
    shell:
        "pathxf {input.config} -i {input.data} "
        "{params.entities}"
