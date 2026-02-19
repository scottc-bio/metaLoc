nextflow.enable.dsl=2

/*
 * =========================
 *        PARAMETERS
 * =========================
 */
params.fasta      = params.fasta
params.organism   = params.organism
params.signalptar = params.signalptar
params.deeploctar = params.deeploctar
params.deeptmhmm = params.deeptmhmm ?: false
params.phobiustar = params.phobiustar
params.meta       = params.meta
params.augustus_model = params.augustus_model
params.min_contig_len = params.min_contig_len

/*
 * =========================
 *     HELPER FUNCTION
 * =========================
 */
def requireFile(path, label) {
    if (!file(path).exists())
        error """
Required ${label} tarball not found: ${path}

This tool requires manual download after license agreement.
See README for instructions.
"""
    return Channel.fromPath(path)
}

/*
 * =========================
 *        PROCESSES
 * =========================
 */

process SETUP_SIGNALP {
    cache 'deep'
    conda 'envs/signalp6.yml'

    input: path tarball

    output: path "assets/signalp6_fast/signalp-6-package/models", emit: models

    script:
    """
    mkdir -p assets
    tar -xzf "$tarball" -C assets
    pip install assets/signalp6_fast/signalp-6-package
    """
}

process RUN_SIGNALP {
    conda 'envs/signalp6.yml'

    publishDir "results/${fasta.baseName}", mode: 'copy'

    input:
        path fasta
        path model_dir
        val organism

    output:
        path "signalp_out/signalp_results.txt", emit: txt
        path "signalp_out/*"

    script:
    """
    mkdir -p signalp_out

    if [ ! -s "$fasta" ]; then
        # Empty but valid SignalP output (no warnings)
        cat > signalp_out/signalp_results.txt <<'EOF'
# SignalP-6.0
# ID	Prediction	CS_Position
EOF
    else
        signalp6 --fastafile "$fasta" --output_dir signalp_out \
                 --organism "$organism" --format txt --mode fast \
                 --model_dir "$model_dir"

        cp signalp_out/prediction_results.txt signalp_out/signalp_results.txt

        tar -czf signalp_out/signalp_plots.tar.gz signalp_out/*_plot.txt
        rm signalp_out/*_plot.txt

        rm signalp_out/prediction_results.txt
    fi
    """
}

process SETUP_DEEPLOC {
    cache 'deep'
    conda 'envs/deeploc2.yml'

    input: path tarball

    output: path "deeploc.installed", emit: installed

    script:
    """
    mkdir -p assets/deeploc2
    tar -xzf "$tarball" -C assets/deeploc2
    pip install assets/deeploc2/deeploc2_package
    touch deeploc.installed
    """
}

process RUN_DEEPLOC {
    conda 'envs/deeploc2.yml'
    publishDir "results/${fasta.baseName}", mode: 'copy'

    input:
        path fasta
        path ready

    output:
        path "deeploc_out/deeploc_results.txt", emit: txt

    script:
    """
    mkdir -p deeploc_out

    if [ ! -s "$fasta" ]; then
        echo "Protein_ID,Localization" > deeploc_out/deeploc_results.txt
    else
        deeploc2 --fasta "$fasta" --output deeploc_out --model Fast
        cp deeploc_out/results_*.csv deeploc_out/deeploc_results.txt
    fi
    """
}

process RUN_DEEPLOCPRO {
    conda 'envs/deeplocpro.yml'

    publishDir "results/${fasta.baseName}", mode: 'copy'

    input: path fasta

    output:
        path "deeplocpro_out/deeplocpro_results.txt", emit: txt
        path "deeplocpro_out/*"

    script:
    """
    mkdir -p deeplocpro_out

    if [ ! -s "$fasta" ]; then
        echo "ACC,Localization" > deeplocpro_out/deeplocpro_results.txt
    else
        sed '/^>/! s/\\*//g' "$fasta" > cleaned.faa
        deeplocpro -f cleaned.faa -o deeplocpro_out
        cp deeplocpro_out/*.csv deeplocpro_out/deeplocpro_results.txt
        mv cleaned.faa deeplocpro_out/cleaned.faa
    fi
    """
}

process RUN_DEEPTMHMM {
    conda 'envs/deeptmhmm.yml'

    publishDir "results/${fasta.baseName}", mode: 'copy'

    input:
        path fasta
        path parser

    output:
        path "deeptmhmm_out/deeptmhmm_results.txt", emit: txt
        path "deeptmhmm_out/*"

    script:
    """
    mkdir -p deeptmhmm_out
    BIO_LIB="\$CONDA_PREFIX/bin/biolib"

    if [ ! -s "$fasta" ]; then
        echo "Sequence_ID,Prediction,Topology" > deeptmhmm_out/deeptmhmm_results.txt
        exit 0
    fi

    sed '/^>/! s/\\*//g' "$fasta" > cleaned.faa
    seqkit seq -m 30 -g cleaned.faa -o cleaned_min30.faa

    "\$BIO_LIB" run DTU/DeepTMHMM --fasta cleaned_min30.faa

    python "$parser" biolib_results/predicted_topologies.3line deeptmhmm_results.txt
    mv deeptmhmm_results.txt deeptmhmm_out/
    mv cleaned_min30.faa deeptmhmm_out/
    mv biolib_results/* deeptmhmm_out/
    """
}

process RUN_PHOBIUS {
    conda 'envs/phobius.yml'
    maxForks 1

    publishDir "results/${fasta.baseName}", mode: 'copy'

    input:
        path fasta
        path tarball

    output:
        path "phobius_out/phobius_results.txt", emit: txt
        path "phobius_out/cleaned.faa", optional: true

    when:
        tarball

    script:
    """

    mkdir -p phobius_out

    sed -E "/^>/! { s/\\*//g; s/[^ACDEFGHIKLMNPQRSTVWYX]/X/g }" "$fasta" > cleaned.faa

    if [ ! -s "$fasta" ]; then
        # Fixed-width header only
        echo "SEQENCE ID TM SP PREDICTION" > phobius_out/phobius_results.txt
        exit 0
    fi

    PHOBIUS_DIR="\$CONDA_PREFIX/share/phobius-1.01-5"

    phobius-register "$tarball"

    PHOBIUS_PL="\$CONDA_PREFIX/share/phobius-1.01-5/phobius.pl"

    perl "\$PHOBIUS_PL" -short cleaned.faa > phobius_out/phobius_results.txt || true

    mv cleaned.faa phobius_out/cleaned.faa
    """
}

process RUN_PRODIGAL {
    conda 'envs/prodigal.yml'

    publishDir "results/prodigal", mode: 'copy'

    input: path contigs

    output:
        path "prok_proteins.faa", emit: proteins
        path "no_contigs_warning.txt", optional: true
        path "no_proteins_warning.txt", optional: true
        path "prok_genes_out.fna", optional: true
        path "prok_genes_out.gff", optional: true

    script:
    """
    if [ ! -s "$contigs" ]; then
        # No contigs at all
        touch prok_proteins.faa
        echo "No prokaryotic contigs provided to Prodigal." > no_contigs_warning.txt
       exit 0
    fi

    prodigal -i "$contigs" -a prok_proteins.faa -d prok_genes_out.fna -f gff -o prok_genes_out.gff -p meta

    if [ ! -s prok_proteins.faa ]; then
        echo "Contigs provided, but no ORFs predicted by Prodigal." \
            > no_proteins_warning.txt
    fi
    """
}

process RUN_AUGUSTUS {
    tag "${augustus_model}"
    conda 'envs/augustus.yml'

    publishDir "results/augustus", mode: 'copy'

    input:
        path contigs
        val augustus_model

    output:
        path "euk_proteins.faa", emit: proteins
        path "no_contigs_warning.txt", optional: true
        path "no_proteins_warning.txt", optional: true
        path "*.fai", optional: true
        path "augustus_genes.gff", optional: true

    script:
    """
    set -euo pipefail

    mkdir -p assets/augustus_config
    cp -r \$CONDA_PREFIX/config/* assets/augustus_config/
    export AUGUSTUS_CONFIG_PATH=\$(realpath assets/augustus_config)

    if [ ! -s "$contigs" ]; then
        # No contigs at all
        touch euk_proteins.faa
        echo "No eukaryotic contigs provided to Augustus." > no_contigs_warning.txt
        exit 0
    fi

    augustus --species="${augustus_model}" \
             --genemodel=partial \
             --strand=both \
             --gff3=on \
             --protein=on \
             --uniqueGeneId=true \
             "$contigs" > augustus_genes.gff

    gffread augustus_genes.gff -g "$contigs" -y euk_proteins.faa

    if [ ! -s euk_proteins.faa ]; then
        echo "Contigs provided, but no protein-coding genes predicted by Augustus." \
            > no_proteins_warning.txt
    fi

    """
}


process RUN_EUKREP {
    tag "min_contig_len = ${min_contig_len}"
    conda 'envs/eukrep.yml'

    publishDir "results/eukrep", mode: 'copy'

    input: 
        path contigs
        val min_contig_len

    output:
        path "euk_contigs.fasta", emit: euk
        path "prok_contigs.fasta", emit: prok
        path "filtered_contigs.fasta", emit: filtered

    script:
    """
    # Filter contigs by length BEFORE classification
    seqkit seq -m ${min_contig_len} "$contigs" > filtered_contigs.fasta
    # Handle edge case: nothing left after filtering
    if [ ! -s filtered_contigs.fasta ]; then
        touch euk_contigs.fasta
        touch prok_contigs.fasta
        exit 0
    fi

    # Run EukRep on length filtered contigs only
    EukRep -i filtered_contigs.fasta \
           -o euk_contigs.fasta \
           --prokarya prok_contigs.fasta
    """
}

process MERGE_RESULTS {
    conda 'envs/merge_results.yml'

    input:
        tuple path(fasta), val(name)
        path txts
        path script

    publishDir "results/${name}_final", mode: 'copy'

    output:
        path "final_merged.tsv"
        file "final_concise.tsv"

    script:
    """
    python "$script" ${txts.join(' ')}
    """
}

/*
 * =========================
 *      PROCESS WRAPPERS
 * =========================
 */

// Wrap AUGUSTUS in a workflow
workflow RUN_AUGUSTUS_WRAPPER {
    take: fasta_ch
    main:
        proteins_ch = RUN_AUGUSTUS(fasta_ch, params.augustus_model).proteins
    emit: proteins_ch
}

// Wrap PRODIGAL in a workflow
workflow RUN_PRODIGAL_WRAPPER {
    take: fasta_ch
    main:
        proteins_ch = RUN_PRODIGAL(fasta_ch).proteins
    emit: proteins_ch
}

/*
 * =========================
 *        SUBWORKFLOW
 * =========================
 */

workflow CORE_PIPELINE {
    take:
        fasta_ch
        organism
        signalp_models_ch
        deeploc_ready_ch
        parser_ch
        phobius_tar_ch
        merge_script_ch
        deeptmhmm_val

    main:
        // Map fasta files to (file, sample_name) tuples
        sample_ch = fasta_ch.map { f -> tuple(f, f.baseName) }

        // --- SignalP always runs ---
        signalp_txt = RUN_SIGNALP(
            fasta_ch,
            signalp_models_ch,
            organism
        ).txt

        // --- DeepLoc / DeepLocPro conditional ---
        def loc_txt
        if (organism == 'euk') {
            loc_txt = RUN_DEEPLOC(fasta_ch, deeploc_ready_ch).txt
        } else {
            loc_txt = RUN_DEEPLOCPRO(fasta_ch).txt
        }

        // --- DeepTMHMM / Phobius conditional ---
        def tm_txt
        if (deeptmhmm_val) {
            // Explicit opt-in
            tm_txt = RUN_DEEPTMHMM(fasta_ch, parser_ch).txt
        } else if (phobius_tar_ch) {
            // Default behavior
            tm_txt = RUN_PHOBIUS(fasta_ch, phobius_tar_ch).txt
        } else {
            tm_txt = Channel.empty()
        }

        // --- Combine all main .txt files ---
        all_txts = signalp_txt
            .combine(loc_txt)
            .combine(tm_txt)

        // --- Merge final results ---
        MERGE_RESULTS(sample_ch, all_txts, merge_script_ch)
}

/*
 * =========================
 *      TOP WORKFLOW
 * =========================
 */

workflow RUN_CORE_EUK {
    take:
        fasta_ch
        signalp_models_ch
        deeploc_installed_ch
        parser_ch
        phobius_tar_ch
        merge_script_ch
        deeptmhmm_val
    main:
        CORE_PIPELINE(
            fasta_ch,
            'euk',
            signalp_models_ch,
            deeploc_installed_ch,
            parser_ch,
            phobius_tar_ch,
            merge_script_ch,
            deeptmhmm_val
        )
}

workflow RUN_CORE_PROK {
    take:
        fasta_ch
        signalp_models_ch
        deeploc_installed_ch
        parser_ch
        phobius_tar_ch
        merge_script_ch
        deeptmhmm_val
    main:
        CORE_PIPELINE(
            fasta_ch,
            'other',
            signalp_models_ch,
            deeploc_installed_ch,
            parser_ch,
            phobius_tar_ch,
            merge_script_ch,
            deeptmhmm_val
        )
}

/*
 * =========================
 *      MAIN WORKFLOW
 * =========================
 */

workflow {

    /*
     * Parameter validation
     */
    if (!params.fasta)
        error "You must supply --fasta"

    if (!params.meta && !params.organism)
        error "--organism is required when --meta is false"

    if ((params.meta in ['euk','mixed']) && !params.augustus_model) {
        error """
        You must supply --augustus_model when --meta is 'euk' or 'mixed'.

        Augustus models are species-specific and strongly affect gene prediction.

        usage:
          --augustus_model=IDENTIFIER

        See the full list of identifiers with:
          cat assets/augustus_model_list.txt
        """
    }

    if (params.meta == 'mixed' && !params.min_contig_len) {
        error """
        When --meta mixed is used, you must supply --min_contig_len.

        EukRep is unreliable on short contigs.
        Values >= 3000 bp are recommended, and values <= 1000 bp are strongly discouraged

        Example:
          --meta mixed --min_contig_len 3000
    """
}

    /*
     * Static inputs (stable cache keys)
     */
    fasta_ch        = Channel.fromPath(params.fasta)
    parser_ch       = Channel.fromPath('bin/deeptmhmm_parse.py')
    merge_script_ch = Channel.fromPath('bin/merge_results.py')

    signalp_tar_ch  = requireFile(params.signalptar, 'SignalP')

    deeploc_tar_ch  = (params.organism == 'euk' || params.meta in ['euk','mixed'])
        ? requireFile(params.deeploctar, 'DeepLoc2')
        : Channel.empty()

    phobius_tar_ch = requireFile(params.phobiustar, 'Phobius')

    /*
     * GLOBAL TOOL SETUP (run once, cacheable)
     */
    signalp_models_ch = SETUP_SIGNALP(signalp_tar_ch).models

    deeploc_installed_ch = (params.organism == 'euk' || params.meta in ['euk','mixed'])
        ? SETUP_DEEPLOC(deeploc_tar_ch).installed
        : Channel.empty()

    /*
     * Routing logic
     */
    if (!params.meta || params.meta == false) {
        CORE_PIPELINE(
            fasta_ch,
            params.organism,
            signalp_models_ch,
            deeploc_installed_ch,
            parser_ch,
            phobius_tar_ch,
            merge_script_ch,
            params.deeptmhmm
        )

    }
    else if (params.meta == 'euk') {

        proteins_ch = RUN_AUGUSTUS(fasta_ch, params.augustus_model).proteins

        CORE_PIPELINE(
            proteins_ch,
            'euk',
            signalp_models_ch,
            deeploc_installed_ch,
            parser_ch,
            phobius_tar_ch,
            merge_script_ch,
            params.deeptmhmm
        )

    }
    else if (params.meta == 'other') {

        proteins_ch = RUN_PRODIGAL(fasta_ch).proteins

        CORE_PIPELINE(
            proteins_ch,
            'other',
            signalp_models_ch,
            deeploc_installed_ch,
            parser_ch,
            phobius_tar_ch,
            merge_script_ch,
            params.deeptmhmm
        )

    }
    else if (params.meta == 'mixed') {
        // Split contigs with EukRep
        eukrep_out = RUN_EUKREP(fasta_ch, params.min_contig_len)

        // Run eukaryotic branch
        euk_proteins_ch = RUN_AUGUSTUS_WRAPPER(eukrep_out.euk)

        // Run prokaryotic branch
        prok_proteins_ch = RUN_PRODIGAL_WRAPPER(eukrep_out.prok)

        RUN_CORE_EUK(euk_proteins_ch, signalp_models_ch, deeploc_installed_ch,
                     parser_ch, phobius_tar_ch, merge_script_ch, params.deeptmhmm)

        RUN_CORE_PROK(prok_proteins_ch, signalp_models_ch, deeploc_installed_ch,
                      parser_ch, phobius_tar_ch, merge_script_ch, params.deeptmhmm)
    }
}