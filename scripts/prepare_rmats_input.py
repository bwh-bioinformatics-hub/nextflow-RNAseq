import os
import sys


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception("Please provide phenotype file!")
    phenotype_name = sys.argv[1]

    sample_dict = {}
    for _, _, files in os.walk("./"):
        for filename in files:
            # STAR results bed filenames are like "SAMPLE_ID.Aligned.sortedByCoord.out.bam"
            if ".bam" in filename:
                sample_id = filename.split(".")[0]
                sample_dict[sample_id] = filename

    # map each phenotype to all related bam files
    phenotype_bam = {}
    for line in open(phenotype_name).readlines():
        group, sample_id = line.strip().split(",")
        group = group.lower()
        if group not in phenotype_bam:
            phenotype_bam[group] = sample_dict[sample_id]
        else:
            phenotype_bam[group] += ","
            phenotype_bam[group] += sample_dict[sample_id]

    # check if "CTRL" in all phenotypes
    # if not, treat first phenotype as "CTRL"
    if "ctrl" in phenotype_bam:
        with open("ctrl_rmats_input.txt", "w") as f:
            f.write(phenotype_bam["ctrl"])
        for each in phenotype_bam:
            if each != "ctrl":
                with open("%s_rmats_input.txt" % each, "w") as f:
                    f.write(phenotype_bam[each])
    else:
        for i, each in enumerate(phenotype_bam):
            if i == 0:
                with open("%s_ctrl_rmats_input.txt" % each, "w") as f:
                    f.write(phenotype_bam[each])
            else:
                with open("%s_rmats_input.txt" % each, "w") as f:
                    f.write(phenotype_bam[each])
