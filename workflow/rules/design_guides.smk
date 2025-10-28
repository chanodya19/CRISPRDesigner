### Design and score gRNAs from scratch using provided input files

import subprocess
import os

## Use original regions directly (no subtraction needed)
rule prepare_regions:
    input:
        regions = config["regions"]
    output:
        regions_minus_predesigned = "results/GuideDesign/newRegions.bed"
    shell:
        "cp {input.regions} {output.regions_minus_predesigned}"

## Find all NGG PAM guides in the selected regions
## Java tool throws an error with SKIP_SCORING=true, but this is expected
rule find_all_guides:
    input:
        regions = "results/GuideDesign/newRegions.bed",
        genome_fasta = config["genome_fasta"],
        off_target_bits = config["off_target_bits"]
    output:
        guides = "results/GuideDesign/allGuides.bed"
    params:
        mem = config["java_memory"]
    run:
        try:
            command = f"java -Xmx{params.mem} -jar workflow/scripts/CRISPRDesigner.jar \
TARGETS={input.regions} OUTPUT_DIR=results/GuideDesign/ \
GENOME_FASTA={input.genome_fasta} \
LENIENT=false \
OFF_TARGETS={input.off_target_bits} \
SKIP_PAIRING=true \
DIVIDE_AND_CONQUER=false \
SKIP_SCORING=true"
            print("Running: " + command)
            proc_output = subprocess.check_output(command, shell=True)
        except subprocess.CalledProcessError as exc:
            if 'filteredGuides.bed' in str(exc.output):
                print("Snakemake caught Java Exception: Ignore this error, which is expected")
                pass
            else:
                raise

## Split all guides into chunks for parallel scoring
checkpoint scatter_scores:
    input:
        guides = "results/GuideDesign/allGuides.bed"
    output:
        scatterdir = directory('results/GuideDesign/scatter')
    params:
        split_guides = config["split_guides"]
    shell:
        """
        mkdir {output}
        split --lines={params.split_guides} -d {input.guides} {output}/guides.
        for file in {output}/guides.*; do 
          mkdir -p {output}/dir-$(basename $file)
          mv $file {output}/dir-$(basename $file)/allGuides.bed
          mv {output}/dir-$(basename $file) $file
        done
        """

## Score guides using Java tool (off-target scoring)
## SKIP_GENERATION=true means guides are already present
rule score_guides:
    input:
        guide_scatter = 'results/GuideDesign/scatter/guides.{i}/',
        regions = config["regions"],
        genome_fasta = config["genome_fasta"],
        off_target_bits = config["off_target_bits"]
    output:
        scored = 'results/GuideDesign/scatter/guides.{i}/filteredGuides.bed'
    params:
        mem = config["java_memory"],
        cwd = os.getcwd()
    shell:
        """
        (
        cd {input.guide_scatter}
        java -Xmx{params.mem} -jar {params.cwd}/workflow/scripts/CRISPRDesigner.jar \
          TARGETS={params.cwd}/{input.regions} \
          OUTPUT_DIR={params.cwd}/{input.guide_scatter} \
          GENOME_FASTA={input.genome_fasta} \
          LENIENT=false \
          OFF_TARGETS={input.off_target_bits} \
          SKIP_PAIRING=true \
          DIVIDE_AND_CONQUER=false \
          SKIP_SCORING=false \
          SKIP_GENERATION=true
        )
        """

## Aggregate scored guide files from scatter step
def aggregate_input(wildcards):
    checkpoint_output = checkpoints.scatter_scores.get().output.scatterdir
    ivals = glob_wildcards(os.path.join(checkpoint_output, 'guides.{i}/allGuides.bed')).i
    return expand('{dir}/guides.{i}/filteredGuides.bed', dir=checkpoint_output, i=ivals)

## Combine and relabel scored guides
rule gather_and_relabel_guide_scores:
    input:
        aggregate_input
    params:
        genome_sizes = config["genome_sizes"],
        regions = config["regions"]
    output:
        combined = 'results/GuideDesign/filteredGuides.bed'
    shell:
        '''
        cat {input} | \
        awk -v OFS=$'\t' '{{ $14="FILLER"; print $0 }}' | \
        bedtools sort -i stdin -faidx {params.genome_sizes} | \
        uniq | \
        bedtools intersect -a stdin -b {params.regions} -wa -wb | cut -f 1-13,18 > {output.combined}
        '''

## Final filtering: remove guides with "TTTT", keep MIT specificity > 50
## Outputs both a design table and IGV-compatible BED file
rule filter_guides:
    input:
        combined_guides = 'results/GuideDesign/filteredGuides.bed',
        genome_sizes = config["genome_sizes"]
    output:
        design_guides = 'results/GuideDesign/designGuides.txt',
        design_guides_igv = 'results/GuideDesign/designGuides.bed'
    shell:
        """
        echo -e "chr\tstart\tend\tlocus\tscore\tstrand\tGuideSequenceWithPAM\tguideSet\tSSC" > {output.design_guides}
        cat {input.combined_guides} | grep -v "TTTT" | awk '{{if ($5 > 50) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $13 "\t" $14 "\t0" }}' >> {output.design_guides}
        sed 1d {output.design_guides} | cut -f 1-6 | uniq > {output.design_guides_igv}
        """
