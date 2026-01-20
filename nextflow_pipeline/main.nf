#!/usr/bin/env nextflow

/*
========================================================================================
    FIBER-SEQ PROCESSING PIPELINE
========================================================================================
    A Nextflow pipeline for processing PacBio HiFi fiber-seq data through:
    - Alignment (pbmm2)
    - Nucleosome calling (fibertools)
    - QC metrics
    - Pileup generation and bigWig conversion
    - Heatmap visualization (deepTools)
    - Chromosome filtering
    - k-mer phasing (optional)
    - FIRE peak calling (optional)
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check mandatory parameters
if (!params.input) {
    error "ERROR: Missing required parameter --input (samplesheet TSV)"
}
if (!params.fasta) {
    error "ERROR: Missing required parameter --fasta (reference genome FASTA)"
}
if (!params.outdir) {
    error "ERROR: Missing required parameter --outdir (output directory)"
}

/*
========================================================================================
    PARSE SAMPLESHEET
========================================================================================
*/

def parseSamplesheet(samplesheet) {
    Channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: false, sep: '\t')
        .map { row ->
            def sample_name = row[0]
            def raw_bam = file(row[1], checkIfExists: true)
            def tissue = row[2]
            def sequencer = row[3]
            return tuple(sample_name, raw_bam, tissue, sequencer)
        }
}

/*
========================================================================================
    PROCESSES - Part 1: Core Processing
========================================================================================
*/

process PBMM2_INDEX {
    tag "index"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/pbmm2:1.14.99--h9ee0642_0'

    input:
    path fasta

    output:
    path "*.mmi", emit: index

    script:
    def preset = params.pbmm2_preset ?: 'HIFI'
    """
    pbmm2 index \\
        --preset ${preset} \\
        -j ${task.cpus} \\
        ${fasta} \\
        ${fasta.baseName}.mmi
    """
}

process PBMM2_ALIGN {
    tag "${sample_name}"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/pbmm2:1.14.99--h9ee0642_0'

    input:
    tuple val(sample_name), path(raw_bam), val(tissue), val(sequencer)
    path reference

    output:
    tuple val(sample_name), path("${sample_name}.aligned.bam"), val(tissue), val(sequencer), emit: aligned_bam
    path "${sample_name}.aligned.bam.bai", emit: bai

    script:
    def preset = params.pbmm2_preset ?: 'HIFI'
    """
    pbmm2 align \\
        --preset ${preset} \\
        --sort \\
        -j ${task.cpus} \\
        --log-level INFO \\
        ${reference} \\
        ${raw_bam} \\
        ${sample_name}.aligned.bam

    # Create BAI index
    samtools index -@ ${task.cpus} ${sample_name}.aligned.bam
    """
}

process SEQUENCING_QC {
    tag "${sample_name}"
    label 'process_low'

    // Note: Requires SMRT Tools - may need custom container or module load approach
    container params.smrt_container ?: null

    input:
    tuple val(sample_name), path(raw_bam), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("${sample_name}.sequencing_qc_report.pdf"), emit: qc_report, optional: true
    path "${sample_name}.consensusreadset.xml", emit: dataset_xml

    script:
    def smrt_root = params.smrt_root ?: ''
    def pbindex_cmd = smrt_root ? "${smrt_root}/smrtcmds/bin/pbindex" : "pbindex"
    def dataset_cmd = smrt_root ? "${smrt_root}/smrtcmds/bin/dataset" : "dataset"
    def runqc_cmd = smrt_root ? "${smrt_root}/smrtcmds/bin/runqc-reports" : "runqc-reports"
    """
    # Index raw BAM if needed
    if [[ ! -f "${raw_bam}.pbi" ]]; then
        ${pbindex_cmd} -j ${task.cpus} ${raw_bam}
    fi

    # Create ConsensusReadSet
    ${dataset_cmd} create \\
        --type ConsensusReadSet \\
        --force \\
        --name "${sample_name}" \\
        ${sample_name}.consensusreadset.xml \\
        ${raw_bam}

    # Generate Sequencing QC Report
    ${runqc_cmd} -b \\
        --pdf-report ${sample_name}.sequencing_qc_report.pdf \\
        ${sample_name}.consensusreadset.xml || echo "WARNING: QC PDF generation failed, continuing..."
    """
}

process FIBERTOOLS_ADD_NUCLEOSOMES {
    tag "${sample_name}"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/fibertools-rs:0.7.1--h3b373d1_0'

    input:
    tuple val(sample_name), path(aligned_bam), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("${sample_name}.ft.bam"), val(tissue), val(sequencer), emit: ft_bam

    script:
    def ft_options = params.ft_nucleosome_options ?: ''
    """
    ft add-nucleosomes \\
        -t ${task.cpus} \\
        ${ft_options} \\
        ${aligned_bam} \\
        ${sample_name}.ft.bam
    """
}

process SAMTOOLS_INDEX {
    tag "${sample_name}"
    label 'process_low'

    container 'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(sample_name), path(bam), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path(bam), path("${bam}.bai"), val(tissue), val(sequencer), emit: indexed_bam
    path "${bam}.pbi", emit: pbi, optional: true

    script:
    """
    samtools index -@ ${task.cpus} ${bam}

    # Try to create PBI index if pbindex is available
    if command -v pbindex &>/dev/null; then
        pbindex ${bam} || true
    fi
    """
}

process FIBERSEQ_QC {
    tag "${sample_name}"
    label 'process_medium'

    // Note: Requires custom container with tcsh and fiberseq-qc suite
    container params.fiberseq_qc_container ?: null

    input:
    tuple val(sample_name), path(ft_bam), path(bai), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("qc_output/*"), emit: qc_results

    script:
    def fiberseq_qc_root = params.fiberseq_qc_root ?: ''
    def runall_qc = fiberseq_qc_root ? "${fiberseq_qc_root}/src/runall-qc.tcsh" : "runall-qc.tcsh"
    """
    mkdir -p qc_output
    cd qc_output
    tcsh ${runall_qc} . ${sample_name} ../${ft_bam}
    """
}

process FIBERTOOLS_PILEUP {
    tag "${sample_name}"
    label 'process_high'

    container 'https://depot.galaxyproject.org/singularity/fibertools-rs:0.7.1--h3b373d1_0'

    input:
    tuple val(sample_name), path(ft_bam), path(bai), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("${sample_name}.*.pileup.bed"), val(tissue), val(sequencer), emit: pileup

    script:
    def region = params.pileup_region ?: ''
    def ml_threshold = params.ml_threshold ? "--ml ${params.ml_threshold}" : ''
    def ft_options = params.ft_pileup_options ?: ''
    def output_name = region ? "${sample_name}.${region.replaceAll(/[:,-]/, '_')}.pileup.bed" : "${sample_name}.genome_wide.pileup.bed"
    def region_arg = region ? "${ft_bam} ${region}" : "${ft_bam}"
    """
    ft pileup \\
        -t ${task.cpus} \\
        -m -c \\
        ${ft_options} \\
        ${ml_threshold} \\
        ${region_arg} \\
        -o ${output_name}
    """
}

process PILEUP_TO_BEDGRAPH {
    tag "${sample_name}"
    label 'process_low'

    container 'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(sample_name), path(pileup), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("*.bedgraph"), val(tissue), val(sequencer), emit: bedgraphs

    script:
    // Column mapping: coverage=4, fire=5, score=6, nuc=7, msp=8, m6a=9, cpg=10
    def marks = params.marks?.split(',') ?: ['m6a', 'nuc', 'cpg']
    def pileup_columns = [
        'm6a': 9,
        'nuc': 7,
        'cpg': 10,
        'fire': 5,
        'msp': 8,
        'coverage': 4
    ]
    def commands = marks.collect { mark ->
        def col = pileup_columns[mark]
        if (col) {
            """
            awk -F'\\t' -v col=${col} \\
                '/^#/ {next} \$4 > 0 {printf "%s\\t%s\\t%s\\t%.6f\\n", \$1, \$2, \$3, \$col/\$4}' \\
                ${pileup} > ${pileup.baseName}.${mark}.bedgraph
            """
        } else {
            "echo 'WARNING: Unknown mark ${mark}, skipping'"
        }
    }.join('\n')
    """
    ${commands}
    """
}

process BEDGRAPH_TO_BIGWIG {
    tag "${sample_name}"
    label 'process_low'

    // Use module instead of container (ucsc-bedgraphtobigwig not available on Galaxy depot)
    // On HPC: module load ucsc-utils
    container null
    beforeScript 'ml ucsc-utils 2>/dev/null || module load ucsc-utils 2>/dev/null || true'

    input:
    tuple val(sample_name), path(bedgraphs), val(tissue), val(sequencer)
    path chrom_sizes

    output:
    tuple val(sample_name), path("*.bw"), val(tissue), val(sequencer), emit: bigwigs

    script:
    """
    for bg in ${bedgraphs}; do
        mark=\$(basename "\$bg" .bedgraph | rev | cut -d. -f1 | rev)
        output_bw=\$(basename "\$bg" .bedgraph).bw
        bedGraphToBigWig "\$bg" ${chrom_sizes} "\$output_bw" 2>/dev/null || \\
            echo "WARNING: bedGraphToBigWig failed for \$bg"
    done
    """
}

process COMPUTE_MATRIX {
    tag "${sample_name}"
    label 'process_medium'

    // deepTools container from Biocontainers
    container 'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0'

    input:
    tuple val(sample_name), path(bigwigs), val(tissue), val(sequencer)
    path regions_bed

    output:
    tuple val(sample_name), path("*.matrix.gz"), path("*.sorted_regions.bed"), val(tissue), val(sequencer), emit: matrix

    script:
    def mode = params.matrix_mode ?: 'reference-point'
    def upstream = params.upstream ?: 2000
    def downstream = params.downstream ?: 2000
    def ref_point = params.reference_point ?: 'center'
    def body_length = params.region_body_length ?: 5000
    def scale_up = params.scale_upstream ?: 3000
    def scale_down = params.scale_downstream ?: 3000
    def regions_basename = regions_bed.baseName

    def mode_args = mode == 'reference-point' ?
        "reference-point --referencePoint ${ref_point} -a ${downstream} -b ${upstream}" :
        "scale-regions --regionBodyLength ${body_length} -a ${scale_down} -b ${scale_up}"
    """
    for bw in ${bigwigs}; do
        mark=\$(basename "\$bw" .bw | rev | cut -d. -f1 | rev)
        matrix_out="${sample_name}.${regions_basename}.\${mark}.matrix.gz"
        sorted_out="${sample_name}.${regions_basename}.\${mark}.sorted_regions.bed"

        computeMatrix ${mode_args} \\
            -S "\$bw" \\
            -R ${regions_bed} \\
            -p ${task.cpus} \\
            --missingDataAsZero \\
            --skipZeros \\
            -o "\$matrix_out" \\
            --outFileSortedRegions "\$sorted_out"
    done
    """
}

process PLOT_HEATMAP {
    tag "${sample_name}"
    label 'process_low'

    // deepTools container from Biocontainers
    container 'https://depot.galaxyproject.org/singularity/deeptools:3.5.5--pyhdfd78af_0'

    input:
    tuple val(sample_name), path(matrices), path(sorted_regions), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("*.heatmap.png"), path("*.heatmap.pdf"), emit: heatmaps

    script:
    def height = params.heatmap_height ?: 15
    def width = params.heatmap_width ?: 4
    def colormap_m6a = params.heatmap_colormap_m6a ?: 'Greens'
    def colormap_nuc = params.heatmap_colormap_nuc ?: 'Blues'
    def colormap_cpg = params.heatmap_colormap_cpg ?: 'Reds'
    """
    for matrix in ${matrices}; do
        mark=\$(basename "\$matrix" .matrix.gz | rev | cut -d. -f1 | rev)
        png_out=\$(basename "\$matrix" .matrix.gz).heatmap.png
        pdf_out=\$(basename "\$matrix" .matrix.gz).heatmap.pdf

        # Set colormap based on mark (using configurable parameters)
        case "\$mark" in
            m6a) colormap="${colormap_m6a}" ;;
            nuc) colormap="${colormap_nuc}" ;;
            cpg) colormap="${colormap_cpg}" ;;
            *) colormap="RdYlBu_r" ;;
        esac

        # Set label based on mark
        case "\$mark" in
            m6a) label="% 6mA" ;;
            nuc) label="% Nucleosome" ;;
            cpg) label="% 5mC" ;;
            *) label="\$mark" ;;
        esac

        plotHeatmap -m "\$matrix" -o "\$png_out" \\
            --colorMap "\$colormap" \\
            --heatmapHeight ${height} --heatmapWidth ${width} \\
            --plotTitle "${sample_name} - \$label" \\
            --legendLocation best

        plotHeatmap -m "\$matrix" -o "\$pdf_out" \\
            --colorMap "\$colormap" \\
            --heatmapHeight ${height} --heatmapWidth ${width} \\
            --plotTitle "${sample_name} - \$label" \\
            --legendLocation best
    done
    """
}

/*
========================================================================================
    PROCESSES - Part 2: Advanced Analysis
========================================================================================
*/

process SAMTOOLS_FILTER_CHROMS {
    tag "${sample_name}"
    label 'process_medium'

    container 'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1'

    input:
    tuple val(sample_name), path(ft_bam), path(bai), val(tissue), val(sequencer)

    output:
    tuple val(sample_name), path("${sample_name}.filtered.bam"), path("${sample_name}.filtered.bam.bai"), val(tissue), val(sequencer), emit: filtered_bam

    script:
    def chromosomes = params.filter_chromosomes ?: ''
    """
    if [[ -n "${chromosomes}" ]]; then
        # Filter reads to specified chromosomes
        samtools view -@ ${task.cpus} -b -h ${ft_bam} ${chromosomes} > ${sample_name}.tmp.bam

        # Clean header - remove @SQ lines for chromosomes not in our list
        samtools view -H ${sample_name}.tmp.bam > full_header.sam
        grep -v "^@SQ" full_header.sam > clean_header.sam
        for chrom in ${chromosomes}; do
            grep -P "\\tSN:\${chrom}\\t" full_header.sam >> clean_header.sam || true
        done

        # Apply clean header
        samtools reheader clean_header.sam ${sample_name}.tmp.bam > ${sample_name}.filtered.bam
        rm -f ${sample_name}.tmp.bam full_header.sam clean_header.sam
    else
        # No filtering - just copy
        cp ${ft_bam} ${sample_name}.filtered.bam
    fi

    # Index
    samtools index -@ ${task.cpus} ${sample_name}.filtered.bam
    """
}

process KMER_PHASING {
    tag "${sample_name}"
    label 'process_very_high'

    input:
    tuple val(sample_name), path(input_bam), path(bai), val(tissue), val(sequencer)
    path fasta
    path pixi_manifest
    path profile_path

    output:
    tuple val(sample_name), path("output/**"), val(tissue), val(sequencer), emit: phasing_results
    tuple val(sample_name), path("output/*.haplotagged.bam"), path("output/*.haplotagged.bam.bai"), val(tissue), val(sequencer), emit: phased_bam, optional: true

    script:
    def resume_flag = params.resume_kmer_phasing ? '--rerun-incomplete' : ''
    def align_flag = params.kmer_align ? 'true' : 'false'
    def parental_mode = params.kmer_parental_mode
    def maternal_bam = params.kmer_maternal_bam ?: ''
    def paternal_bam = params.kmer_paternal_bam ?: ''
    """
    mkdir -p output
    cd output

    # Create sample-specific config YAML dynamically
    cat > config.yaml << 'CONFIGEOF'
sample: ${sample_name}
bam: ../${input_bam}
ref: ../${fasta}
align: ${align_flag}
CONFIGEOF

    # Add parental data if enabled
    if [[ "${parental_mode}" == "true" ]]; then
        if [[ -n "${maternal_bam}" ]]; then
            echo "maternal: ${maternal_bam}" >> config.yaml
        fi
        if [[ -n "${paternal_bam}" ]]; then
            echo "paternal: ${paternal_bam}" >> config.yaml
        fi
    fi

    # Run snakemake via pixi
    pixi run --manifest-path ../${pixi_manifest} snakemake \\
        --profile ../${profile_path} \\
        --configfile config.yaml \\
        --verbose \\
        ${resume_flag}
    """
}

process FIRE {
    tag "${sample_name}"
    label 'process_very_high'

    input:
    tuple val(sample_name), path(input_bam), path(bai), val(tissue), val(sequencer)
    path fasta
    path pixi_manifest
    path profile_path

    output:
    tuple val(sample_name), path("output/**"), emit: fire_results

    script:
    // Required params
    def ref_name = params.fire_ref_name ?: 'hg38'
    def min_contig_length = params.fire_min_contig_length ?: 0

    // Optional params (null means use FIRE defaults)
    def ont_flag = params.fire_ont ? 'true' : 'false'
    def max_threads = params.fire_max_threads
    def keep_chromosomes = params.fire_keep_chromosomes
    def excludes = params.fire_excludes
    def max_peak_fdr = params.fire_max_peak_fdr
    def min_fire_fdr = params.fire_min_fire_fdr
    def min_coverage = params.fire_min_coverage
    def coverage_within_n_sd = params.fire_coverage_within_n_sd
    def min_msp = params.fire_min_msp
    def min_ave_msp_size = params.fire_min_ave_msp_size
    def min_per_acc_peak = params.fire_min_per_acc_peak
    def min_frac_accessible = params.fire_min_frac_accessible
    """
    mkdir -p output
    cd output

    # Create manifest TSV for FIRE (sample, bam, fai columns)
    echo -e "sample\\tbam\\tfai" > manifest.tsv
    echo -e "${sample_name}\\t../${input_bam}\\t../${bai}" >> manifest.tsv

    # Generate FIRE config YAML dynamically
    cat > fire_config.yaml << 'CONFIGEOF'
# Auto-generated FIRE configuration
ref: ../${fasta}
ref_name: ${ref_name}
manifest: manifest.tsv
min_contig_length: ${min_contig_length}
ont: ${ont_flag}
CONFIGEOF

    # Add optional parameters only if they are set (non-null)
    ${max_threads ? "echo 'max_t: ${max_threads}' >> fire_config.yaml" : ''}
    ${keep_chromosomes ? "echo 'keep_chromosomes: \"${keep_chromosomes}\"' >> fire_config.yaml" : ''}
    ${max_peak_fdr ? "echo 'max_peak_fdr: ${max_peak_fdr}' >> fire_config.yaml" : ''}
    ${min_fire_fdr ? "echo 'min_fire_fdr: ${min_fire_fdr}' >> fire_config.yaml" : ''}
    ${min_coverage ? "echo 'min_coverage: ${min_coverage}' >> fire_config.yaml" : ''}
    ${coverage_within_n_sd ? "echo 'coverage_within_n_sd: ${coverage_within_n_sd}' >> fire_config.yaml" : ''}
    ${min_msp ? "echo 'min_msp: ${min_msp}' >> fire_config.yaml" : ''}
    ${min_ave_msp_size ? "echo 'min_ave_msp_size: ${min_ave_msp_size}' >> fire_config.yaml" : ''}
    ${min_per_acc_peak ? "echo 'min_per_acc_peak: ${min_per_acc_peak}' >> fire_config.yaml" : ''}
    ${min_frac_accessible ? "echo 'min_frac_accessible: ${min_frac_accessible}' >> fire_config.yaml" : ''}

    # Handle excludes list if provided
    ${excludes ? """
    echo 'excludes:' >> fire_config.yaml
    for excl in ${excludes}; do
        echo "  - \\\"\$excl\\\"" >> fire_config.yaml
    done
    """ : ''}

    # Run FIRE via pixi
    pixi run --manifest-path ../${pixi_manifest} fire \\
        --profile ../${profile_path} \\
        --configfile fire_config.yaml
    """
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {

    // Parse samplesheet
    ch_samples = parseSamplesheet(params.input)

    // Reference handling
    ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)

    // Optionally create or use existing index
    if (params.fasta_index) {
        ch_reference = Channel.fromPath(params.fasta_index, checkIfExists: true)
    } else {
        PBMM2_INDEX(ch_fasta)
        ch_reference = PBMM2_INDEX.out.index
    }

    /*
    ----------------------------------------------------------------------------------------
        Part 1: Core Processing
    ----------------------------------------------------------------------------------------
    */

    // Step 1: Alignment (skip if pre-aligned BAMs provided)
    if (params.skip_alignment) {
        // bam_path in samplesheet is already aligned - pass directly
        log.info "Skipping alignment step - using pre-aligned BAMs from samplesheet"
        if (!params.fasta_index) {
            log.warn "skip_alignment is set but no fasta_index provided - reference indexing will still run (unnecessary)"
        }
        ch_aligned = ch_samples
    } else {
        PBMM2_ALIGN(ch_samples, ch_reference.collect())
        ch_aligned = PBMM2_ALIGN.out.aligned_bam
    }

    // Step 1b: Sequencing QC (runs on raw BAM, parallel with alignment - only for raw BAMs)
    if (params.run_sequencing_qc && !params.skip_alignment) {
        SEQUENCING_QC(ch_samples)
    }

    // Step 2: Add nucleosomes
    FIBERTOOLS_ADD_NUCLEOSOMES(ch_aligned)

    // Step 2b: Index ft.bam
    SAMTOOLS_INDEX(FIBERTOOLS_ADD_NUCLEOSOMES.out.ft_bam)

    // Step 3: Fiber-seq QC (optional, requires special container)
    if (params.run_fiberseq_qc) {
        FIBERSEQ_QC(SAMTOOLS_INDEX.out.indexed_bam)
    }

    // Step 4: Create pileups
    FIBERTOOLS_PILEUP(SAMTOOLS_INDEX.out.indexed_bam)

    // Step 5a: Convert pileup to bedGraph
    PILEUP_TO_BEDGRAPH(FIBERTOOLS_PILEUP.out.pileup)

    // Step 5b: Convert bedGraph to bigWig
    if (params.chrom_sizes) {
        ch_chrom_sizes = Channel.fromPath(params.chrom_sizes, checkIfExists: true)
        BEDGRAPH_TO_BIGWIG(PILEUP_TO_BEDGRAPH.out.bedgraphs, ch_chrom_sizes.collect())

        // Step 6: Create heatmaps (only if regions_bed provided)
        if (params.regions_bed) {
            ch_regions = Channel.fromPath(params.regions_bed, checkIfExists: true)
            COMPUTE_MATRIX(BEDGRAPH_TO_BIGWIG.out.bigwigs, ch_regions.collect())
            PLOT_HEATMAP(COMPUTE_MATRIX.out.matrix)
        }
    }

    /*
    ----------------------------------------------------------------------------------------
        Part 2: Advanced Analysis (Optional)
    ----------------------------------------------------------------------------------------
    */

    // Step 7: Filter chromosomes (if filtering is requested or needed for downstream)
    if (params.filter_chromosomes || params.run_kmer_phasing || params.run_fire) {
        SAMTOOLS_FILTER_CHROMS(SAMTOOLS_INDEX.out.indexed_bam)
    }

    // Step 8: k-mer phasing (optional)
    if (params.run_kmer_phasing && params.kmer_pixi_manifest) {
        ch_pixi_kmer = Channel.fromPath(params.kmer_pixi_manifest, checkIfExists: true)
        ch_profile_kmer = params.kmer_snakemake_profile ?
            Channel.fromPath(params.kmer_snakemake_profile, checkIfExists: true) :
            Channel.empty()

        // Choose input BAM based on kmer_use_filtered parameter
        ch_kmer_input = params.kmer_use_filtered && params.filter_chromosomes ?
            SAMTOOLS_FILTER_CHROMS.out.filtered_bam :
            SAMTOOLS_INDEX.out.indexed_bam

        KMER_PHASING(
            ch_kmer_input,
            ch_fasta.collect(),
            ch_pixi_kmer.collect(),
            ch_profile_kmer.collect()
        )
    }

    // Step 9: FIRE peak calling (optional)
    if (params.run_fire && params.fire_pixi_manifest) {
        ch_pixi_fire = Channel.fromPath(params.fire_pixi_manifest, checkIfExists: true)
        ch_profile_fire = params.fire_snakemake_profile ?
            Channel.fromPath(params.fire_snakemake_profile, checkIfExists: true) :
            Channel.empty()

        // Choose input BAM based on fire_input_source parameter
        // Options: "nucleosomes" (step 2/3) or "kmer_phasing" (step 8)
        if (params.fire_input_source == 'kmer_phasing' && params.run_kmer_phasing) {
            // Use phased BAM from k-mer phasing
            ch_fire_input = KMER_PHASING.out.phased_bam
        } else if (params.filter_chromosomes) {
            // Use filtered BAM from step 7
            ch_fire_input = SAMTOOLS_FILTER_CHROMS.out.filtered_bam
        } else {
            // Use nucleosome-called BAM from step 2/3
            ch_fire_input = SAMTOOLS_INDEX.out.indexed_bam
        }

        FIRE(
            ch_fire_input,
            ch_fasta.collect(),
            ch_pixi_fire.collect(),
            ch_profile_fire.collect()
        )
    }
}

/*
========================================================================================
    WORKFLOW INTROSPECTION
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed at: ${workflow.complete}"
    log.info "Execution status: ${workflow.success ? 'OK' : 'FAILED'}"
    log.info "Duration: ${workflow.duration}"
    log.info ""
}

workflow.onError {
    log.error "Pipeline execution stopped with error: ${workflow.errorMessage}"
}
