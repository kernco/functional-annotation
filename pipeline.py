import subprocess
import sys
import os
import shlex
import json
import re


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
        self.process = subprocess.Popen(shlex.split(self.command), stdin=stdin, stdout=stdout)
        return self.process


class Pipeline:
    def __init__(self, jsonstr):
        pipeline = json.loads(jsonstr)
        self.name = pipeline["name"]
        self.commands = [task["command"] for task in pipeline["tasks"]]
        self.tasks = []

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

    def setup(self):
        parameters = set()
        for command in self.commands:
            parameters.update(re.findall(r"{(.*?)}", command))
        with open("config.txt", "w") as outfile:
            outfile.write("{\n")
            paramlist = list(parameters)
            for parameter in paramlist[:-1]:
                outfile.write('\t"{}": "",\n'.format(parameter))
            outfile.write('\t"{}": ""\n'.format(paramlist[-1]))
            outfile.write("}\n")

    def load_parameters(self, configfile):
        with open(configfile) as f:
            self.parameters = json.loads(f.read())

    def initialize(self):
        self.tasks = []
        for command in self.commands:
            command = command.format(**self.parameters)
            commands = command.split('|')
            command, infile, outfile = self._parse_command(commands[0])
            self.tasks.append(Task(command, infile=infile, outfile=outfile))
            if len(commands) > 1:
                self.tasks[-1].outfile = subprocess.PIPE
                for command, infile, outfile in [self._parse_command(x) for x in commands[1:-1]]:
                    self.tasks.append(Task(command, infile=self.tasks[-1], outfile=subprocess.PIPE))
                command, infile, outfile = self._parse_command(commands[-1])
                self.tasks.append(Task(command, infile=self.tasks[-1], outfile=outfile))

    def dry_run(self):
        for task in self.tasks:
            print task

    def run(self):
        if not self.tasks:
            self.initialize()
        for task in self.tasks:
            proc = task.run()
            if proc.stdout != subprocess.PIPE:
                proc.wait()


if __name__ == "__main__":
    with open(sys.argv[2]) as f:
        pipeline = Pipeline(f.read())
    if sys.argv[1] == "Setup":
        pipeline.setup()
        print "Pipeline configuration created"
        print "Edit config.txt to specify pipeline parameters"
    elif sys.argv[1] == "Test":
        pipeline.load_parameters("config.txt")
        pipeline.initialize()
        pipeline.dry_run()
    elif sys.argv[1] == "Run":
        pipeline.load_parameters("config.txt")
        pipeline.run()

