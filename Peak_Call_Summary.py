import sys
import subprocess

assay = sys.argv[1]
tissues = sys.argv[2].split()
reps = sys.argv[3].split()

print "{: >12s}".format("Tissue") + ''.join(["{: >12s}".format(rep) for rep in reps]) + "{: >12s}".format("Combined")
for tissue in tissues:
    outline = "{: >12s}".format(tissue)
    for rep in reps:
        try:
            output = subprocess.check_output(['wc', '-l', 'Peak_Calls/{}_{}_{}/macs2_peaks.bed'.format(assay, tissue, rep)])
            outline += "{: >12,d}".format(int(output.split()[0]))
        except subprocess.CalledProcessError:
            outline += "{: >12s}".format("N/A")
    try:
        output = subprocess.check_output(['wc', '-l', 'Peak_Calls/{}_{}_Combined_Peaks.bed'.format(assay, tissue)])
        outline += "{: >12,d}".format(int(output.split()[0]))
    except subprocess.CalledProcessError:
        outline += "{: >12s}".format("N/A")
    print outline

