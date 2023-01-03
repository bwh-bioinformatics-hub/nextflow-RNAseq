import subprocess
import glob
import sys


def getstatusoutput(cmd):
    """ behave like commands.getstatusoutput which moved to
    subprocess.getstatusoutput in Python 3.
    Implmented with subprocess.check_output which is available
    in both Python 2 and 3.
    """
    status = 0
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        status = e.returncode
        output = e.output

    output = output.decode()
    if output[-1] == '\n':
        output = output[:-1]

    return status, output


if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise Exception("Parameters not enough! Run as: python3 execute_rmats.py GTF_PATH N_THREADS")

    gtf = sys.argv[1]
    n_threads = sys.argv[2]

    files = glob.glob("*_rmats_input.txt")
    if not files:
        raise Exception("input file for rMATs not found!")

    ctrl = ""
    ctrl_sample_name = ""
    for filename in files:
        if "ctrl" in filename:
            ctrl = filename
            ctrl_sample_name = filename.split("_")[0]
    if ctrl == "":
        raise Exception("CTRL Group not found!")

    cmd = "/data/bioinformatics/tools/run_rmats " \
            "-t paired --variable-read-length --allow-clipping --nthread 4 " \
            "--gtf " + gtf + \
            " --b1 " + ctrl
    for filename in files:
        if "ctrl" not in filename:
            sample_name = filename.split("_")[0]
            cur_cmd = cmd + \
                      " --b2 " + filename +  \
                      " --od ./rmats-results-%s-vs-%s" % (sample_name, ctrl_sample_name) + \
                      " --tmp ./tmp-%s-vs-%s" % (sample_name, ctrl_sample_name)

            print(cur_cmd)
            print("RMATS Analyzing: %s vs %s..." % (sample_name, ctrl_sample_name))
            status, output = getstatusoutput(cur_cmd)
            if int(status) != 0:
                print("\n*********** error in analyzing %s vs %s..." % (sample_name, ctrl_sample_name))
                print("exit status: ", status)
                print("error detail: \n%s\n" % output)
                # raise Exception()
            else:
                print(output)
