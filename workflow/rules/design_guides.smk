### Design and score gRNAs

import subprocess
import os

rule find_all_guides:
    ## Find all NGG PAM guides in your enhancer regions
    input:
        regions = config["regions"],
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
        # the Java script will throw an error when using SKIP_SCORING=true, even when everything works


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
        root = config["project_root"]
    shell:
        """
        (
        cd {input.guide_scatter}
        java -Xmx{params.mem} -jar {params.root}/workflow/scripts/CRISPRDesigner.jar \
          TARGETS={params.root}/{input.regions} \
          OUTPUT_DIR={params.root}/{input.guide_scatter} \
          GENOME_FASTA={params.root}/{input.genome_fasta} \
          LENIENT=false \
          OFF_TARGETS={params.root}/{input.off_target_bits} \
          SKIP_PAIRING=true \
          DIVIDE_AND_CONQUER=false \
          SKIP_SCORING=false \
          SKIP_GENERATION=true
        )
        """
        # SKIP_GENERATION=true means that the script will read a set of guides present in "OUTPUT_DIR/allGuides.bed"
        # SKIP_SCORING=false means that the script will conduct the off-targeting scoring calculation


def aggregate_input(wildcards):
    '''
    aggregate the file names of the files generated at the scatter step
    '''
    checkpoint_output = checkpoints.scatter_scores.get().output.scatterdir
    ivals = glob_wildcards(os.path.join(checkpoint_output, 'guides.{i}/allGuides.bed')).i
    #print("ivals={}".format(ivals))
    return expand('{dir}/guides.{i}/filteredGuides.bed', dir=checkpoint_output, i=ivals)


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

rule rename_filtered_guides:
    input:
        newguides = 'results/GuideDesign/filteredGuides.new.bed'
    output:
        combinedGuides = 'results/GuideDesign/filteredGuides.bed'
    shell:
        "cp {input.newguides} {output.combinedGuides}"


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
 
