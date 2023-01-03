import os


def trans(rna_path="./"):
    for _, dirs, _ in os.walk(rna_path):
        for dir_name in dirs:
            if "CTRL" not in dir_name:
                continue
            temp_path = rna_path + dir_name
            for _, _, filenames in os.walk(temp_path):
                for f in filenames:
                    if "whole_differential_expression.txt" in f:
                        file_path = temp_path + "/" + f
                        rna_file = open(file_path)
                        csv_name_1 = temp_path + "/de_up.csv"
                        csv_name_2 = temp_path + "/de_down.csv"

                        csv_file_1 = open(csv_name_1, "w")
                        csv_file_2 = open(csv_name_2, "w")

                        csv_file_1.write("gene_id,log2FC,p_value\n")
                        csv_file_2.write("gene_id,log2FC,p_value\n")

                        for line in rna_file.readlines()[1:]:
                            info = line.strip().split()
                            rna_id = info[0]
                            log2_fold_change = float(info[2])

                            if info[-2] == "NA":
                                p_value = 1.0
                            else:
                                p_value = str(info[-2])

                            value = "%s,%f,%s\n" % (rna_id, log2_fold_change, p_value)
                            if log2_fold_change > 0:
                                csv_file_1.write(value)
                            else:
                                csv_file_2.write(value)

                        csv_file_1.close()
                        csv_file_2.close()

    print("Trans DE results to csv done!")


if __name__ == '__main__':
    trans()
