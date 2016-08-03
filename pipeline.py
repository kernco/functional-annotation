import ruffus
import subprocess




starting_files = ['Raw_Reads/H3K4me3_Liver_A.fq.gz', 'Raw_Reads/H3K4me3_Liver_B.fq.gz']

pipeline = ruffus.Pipeline(name = "test")
pipeline.transform(task_func  = trim_reads,
                   input      = starting_files,
                   filter     = ruffus.suffix(".fq.gz"),
                   output     = ".trimmed.fq.gz",
                   output_dir = "Trimmed_Reads")


parser = ruffus.cmdline.get_argparse(description="")
parser.add_argument("--input_file")
options = parser.parse_args()
logger, logger_mutex = ruffus.cmdline.setup_logging(__name__, options.log_file, options.verbose)


def trim_reads(input_file, output_file):
    #subprocess.call(['trim_galore', '--paired', input_files])
    print "Will trim", input_file, " and produce", output_file

ruffus.cmdline.run(options)
#pipeline.run()
