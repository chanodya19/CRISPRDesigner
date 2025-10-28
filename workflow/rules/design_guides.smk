### Design and score gRNAs

import subprocess
import os

## Combine predesigned guide regions into one BED file
## This rule is safe even if no predesigned guides are provided
rule combine_predesigned_guides:
    input:
        regions = lambda wildcards: [
            config["predesigned_guides"][guide_set]["regions"]
            for guide_set in config.get("predesigned_guides", {})
        ],
        genome_sizes = config["genome_sizes"]
    output:
        combined_regions = "results/GuideDesign/predesignedGuideRegions.bed"
    run:
        if len(input.regions) > 0:
            shell(f"cat {' '.join(input.regions)} | bedtools sort -i stdin -faidx {input.genome_sizes} > {output.combined_regions}")
        else:
            shell(f"touch {output.combined_regions}")

## Filter predesigned guides to those overlapping target regions
## Adds dummy SSC column and intersects with regions
rule get_predesigned_guides:
    input:
        predesigned_guides = lambda wildcards: [
            config["predesigned_guides"][guide_set]["guides"]
            for guide_set in config.get("predesigned_guides", {})
        ],
        regions = config["regions"],
        genome_sizes = config["genome_sizes"]
    output:
        guides_in_regions = "results/GuideDesign/filteredGuides.predesigned.bed"
    run:
        if len(input.predesigned_guides) > 0:
            shell("""
                cat {input.predesigned_guides} | \
                awk -v OFS=$'\t' '{{ $14="FILLER"; print $0 }}' | \
                bedtools sort -i stdin -faidx {input.genome_sizes} | uniq | \
                bedtools intersect -a stdin -b {input.regions} -wa -wb | cut -f 1-13,18 > {output.guides_in_regions}
            """)
        else:
            shell(f"touch {output.guides_in_regions}")

## Subtract predesigned regions from target regions and merge remaining
rule subtract_predesigned_regions_and_merge:
    input:
        predesigned = "results/GuideDesign/predesignedGuideRegions.bed",
        regions = config["regions"],
        genome_sizes = config["genome_sizes"]
    output:
        regions_minus_predesigned = "results/GuideDesign/newRegions.bed"
    shell:
        "bedtools slop -i {input.predesigned} -b -22 -g {input.genome_sizes} |"
        "bedtools sort -i stdin -faidx {input.genome_sizes} |"
        "bedtools subtract -a {input.regions} -b stdin -g {input.genome_sizes} |"
        "bedtools sort -i stdin -faidx {input.genome_sizes} |"
        "bedtools merge -i stdin -g {input.genome_sizes} >"
        "{output.regions_minus_predesigned}"

## Find all NGG PAM guides in the selected regions
## Java tool will throw an error with SKIP_SCORING=true, but this is expected
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
## Adds dummy SSC column and intersects with regions
rule gather_and_relabel_guide_scores:
    input:
        aggregate_input
    params:
        genome_sizes = config["genome_sizes"],
        regions = config["regions"]
    output:
        combined = 'results/GuideDesign/filteredGuides.new.bed'
    shell:
        '''
        cat {input} | \
        awk -v OFS=$'\t' '{{ $14="FILLER"; print $0 }}' | \
        bedtools sort -i stdin -faidx {params.genome_sizes} | \
        uniq | \
        bedtools intersect -a stdin -b {params.regions} -wa -wb | cut -f 1-13,18 > {output.combined}
        '''

## Combine new and predesigned guides (if any)
## If predesigned is empty, this still works safely
rule combine_new_predesigned:
    input:
        newguides = 'results/GuideDesign/filteredGuides.new.bed',
        predesigned = 'results/GuideDesign/filteredGuides.predesigned.bed',
        genome_sizes = config["genome_sizes"]
    output:
        combinedGuides = 'results/GuideDesign/filteredGuides.bed'
    shell:
        """
        cat {input.newguides} {input.predesigned} | bedtools sort -i stdin -faidx {input.genome_sizes} | uniq > {output.combinedGuides}
        """

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
        sed 1d {output.design_guides} |
