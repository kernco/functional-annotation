import subprocess
import sys
import os
import shlex
import json


class Task:
    def __init__(self, command, infile=None, outfile=None, dependencies=None, products=None):
        self.command = command
        self.infile = infile
        self.outfile = outfile
        self.dependencies = dependencies
        self.products = products

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
        if not self.dependencies or not self.products:
            return True # Always run if files not specified
        intime = max([os.path.getmtime(filename) for filename in self.dependencies])
        try:
            outtime = min([os.path.getmtime(filename) for filename in self.products])
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
    def __init__(self, jsonstr):
        pipeline = json.loads(jsonstr)
        self.name = pipeline["name"]
        self.commands = [task["command"] for task in pipeline["tasks"]]

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


with open(sys.argv[1]) as f:
    map_reads = Pipeline(f.read())
map_reads.run({'sample': sys.argv[2], 'threads': 2, 'assembly': sys.argv[3]})

