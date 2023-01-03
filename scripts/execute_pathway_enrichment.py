import subprocess
import os

if __name__ == '__main__':
    for _, _, filenames in os.walk("./results"):
        for filename in filenames:
            if "csv" in filename:
                cmd = "Rscript pathway_enrich.R ./results/" + filename
                print(cmd)
                status, output = subprocess.getstatusoutput(cmd)
                if status != 0:
                    print("\n\nError when enrichment for:", filename)
                    print("Error Message:")
                    print(output)
                    print("\n\n")
