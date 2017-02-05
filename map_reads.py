import subprocess
import sys
import os
import shlex


class Task:
    def __init__(self, command, infile=None, outfile=None, infiles=None, outfiles=None):
        self.command = command
        self.infile = infile
        self.outfile = outfile
        self.infiles = infiles
        self.outfiles = outfiles

    def __str__(self):
        rstr = self.command
        if isinstance(self.infile, Task):
            rstr = '| ' + rstr
        elif self.infile:
            rstr = rstr + " < " + self.infile
        if self.outfile == subprocess.PIPE:
            rstr = rstr + ' |'
        elif self.outfile:
            rstr = rstr + " > " + self.outfile
        return rstr

    def outdated(self):
        if not self.infiles or not self.outfiles:
            return True # Always run if files not specified
        intime = max([os.path.getmtime(filename) for filename in self.infiles])
        try:
            outtime = min([os.path.getmtime(filename) for filename in self.outfiles])
        except OSError:
            return True #An output file doesn't exist
        return intime > outtime #At least one input file is newer than at least one output file

    def run(self):
        stdin = None
        stdout = self.outfile
        if isinstance(self.infile, Task):
            stdin = self.infile.process.stdout
        elif self.infile:
            stdin = open(self.infile, 'r')
        if self.outfile and self.outfile != subprocess.PIPE:
            stdout = open(self.outfile, 'wb')
        #print self.command.format(**params), stdin, stdout

        self.process = subprocess.Popen(shlex.split(self.command), stdin=stdin, stdout=stdout)
        return self.process


class Pipeline:
    def __init__(self, name):
        self.name = name
        self.commands = []

    def add_step(self, command):
        self.commands.append(command)

    def _parse_command(self, command):
        infile = ''
        outfile = ''
        if '<' in command:
            command, infile = command.split('<')
        if '>' in command:
            command, outfile = command.split('>')
        if '<' in outfile:
            outfile, infile = outfile.split('<')
        if '>' in infile:
            infile, outfile = infile.split('>')
        infile = infile.strip()
        outfile = outfile.strip()
        if infile == '':
            infile = None
        if outfile == '':
            outfile = None
        return command.strip(), infile, outfile

    def run(self, params):
        for command in self.commands:
            tasks = []
            command = command.format(**params)
            print command
            commands = command.split('|')
            command, infile, outfile = self._parse_command(commands[0])
            tasks.append(Task(command, infile=infile, outfile=outfile))
            if len(commands) > 1:
                tasks[0].outfile = subprocess.PIPE
                for command, infile, outfile in [self._parse_command(x) for x in commands[1:-1]]:
                    tasks.append(Task(command, infile=tasks[-1], outfile=subprocess.PIPE))
                command, infile, outfile = self._parse_command(commands[-1])
                tasks.append(Task(command, infile=tasks[-1], outfile=outfile))
            for task in tasks:
                print task
                proc = task.run()
            proc.wait()


map_reads = Pipeline("Map Reads")
map_reads.add_step("trim_galore Samples/{sample}/raw_reads.fq.gz -o Samples/{sample}")
map_reads.add_step("bwa aln -t {threads} {assembly} Samples/{sample}/raw_reads_trimmed.fq.gz > Samples/{sample}/bwa_aln.sai")
map_reads.add_step("bwa samse {assembly} Samples/{sample}/bwa_aln.sai Samples/{sample}/raw_reads_trimmed.fq.gz | samtools view -bS - > Samples/{sample}/bwa_aln.bam")
map_reads.add_step("samtools view -F 1804 -q 30 -b Samples/{sample}/bwa_aln.bam > Samples/{sample}/filter1.bam")
map_reads.add_step("samtools sort -@ {threads} -o Samples/{sample}/sorted.bam Samples/{sample}/filter1.bam")
map_reads.add_step("picard-tools MarkDuplicates INPUT=Samples/{sample}/sorted.bam OUTPUT=Samples/{sample}/dup_marked.bam METRICS_FILE=Samples/{sample}/dup_metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false")
map_reads.add_step("samtools view -F 1804 -q 30 -b Samples/{sample}/dup_marked.bam > Samples/{sample}/filter2.bam")
map_reads.run({'sample': sys.argv[1], 'threads': 2, 'assembly': sys.argv[2]})

