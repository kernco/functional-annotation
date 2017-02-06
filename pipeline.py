import subprocess
import sys
import os
import shlex
import json
import re
import logging

class DependencyError(Exception):
    def __init__(self, missing):
        self.missing = missing

class Task:
    def __init__(self, command, name=None, infile=None, outfile=None, dependencies=None, products=None, workdir=None):
        self.command = command
        self.name = name
        self.infile = infile
        self.outfile = outfile
        self.dependencies = dependencies
        self.products = products
        self.workdir = workdir

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
        intime = max([os.path.getmtime(os.path.join(self.workdir, filename)) for filename in self.dependencies])
        try:
            outtime = min([os.path.getmtime(os.path.join(self.workdir, filename)) for filename in self.products])
        except OSError:
            return True #An output file doesn't exist
        return intime > outtime #At least one input file is newer than at least one output file

    def check_dependencies(self):
        if self.dependencies:
            missing = []
            for filename in self.dependencies:
                if not os.path.isfile(os.path.join(self.workdir, filename)):
                    missing.append(filename)
            if missing:
                raise DependencyError(missing)

    def run(self):
        self.check_dependencies()
        if not self.outdated():
            return None
        stdin = None
        stdout = self.outfile
        if isinstance(self.infile, Task):
            stdin = self.infile.process.stdout
        elif self.infile:
            stdin = open(os.path.join(self.workdir, self.infile), 'r')
        if self.outfile and self.outfile != subprocess.PIPE:
            stdout = open(os.path.join(self.workdir, self.outfile), 'wb')
        self.process = subprocess.Popen(shlex.split(self.command), stdin=stdin, stdout=stdout, cwd=self.workdir)
        return self.process


class Pipeline:
    def __init__(self, jsonstr):
        self.parameters = json.loads(jsonstr)
        logging.info("Loading pipeline from {}".format(self.parameters["pipeline"]))
        with open(self.parameters["pipeline"]) as f:
            pipeline = json.loads(f.read())
        self.name = pipeline["name"]
        self.steps = pipeline["tasks"]
        self.workdir = pipeline["workdir"].format(**self.parameters)
        self._initialize()
        self.check_dependencies()

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

    def _logmessage(self, message):
        logging.info("{} in {}: {}".format(self.name, self.workdir, message))

    def parameters(self):
        parameters = set()
        for command in self.commands:
            parameters.update(re.findall(r"{(.*?)}", command))
        return list(parameters)

    def _initialize(self):
        self._logmessage("Initializing")
        self.tasks = []
        for step in self.steps:
            commandline = step["command"].format(**self.parameters)
            if "depends" not in step:
                dependencies = None
            elif not isinstance(step["depends"], list):
                dependencies = [step["depends"].format(**self.parameters)]
            else:
                dependencies = [dependency.format(**self.parameters) for dependency in step["depends"]]
            if "produces" not in step:
                products = None
            elif not isinstance(step["produces"], list):
                products = [step["produces"].format(**self.parameters)]
            else:
                products = [product.format(**self.parameters) for product in step["produces"]]
            commands = commandline.split('|')
            command, infile, outfile = self._parse_command(commands[0])
            self.tasks.append(Task(command, infile=infile, outfile=outfile, dependencies=dependencies, products=products, workdir=self.workdir))
            if len(commands) > 1:
                self.tasks[-1].outfile = subprocess.PIPE
                for command, infile, outfile in [self._parse_command(x) for x in commands[1:-1]]:
                    self.tasks.append(Task(command, infile=self.tasks[-1], outfile=subprocess.PIPE, dependencies=dependencies, products=products, workdir=self.workdir))
                command, infile, outfile = self._parse_command(commands[-1])
                self.tasks.append(Task(command, name=step["name"], infile=self.tasks[-1], outfile=outfile, dependencies=dependencies, products=products, workdir=self.workdir))
            else:
                self.tasks[-1].name = step["name"]

    #Return a list of all the filenames listed as dependencies of a task but
    #not as the product of another task.
    def dependencies(self):
        task_dependencies = set()
        task_products = set()
        for task in self.tasks:
            if task.dependencies:
                task_dependencies.update(task.dependencies)
            if task.products:
                task_products.update(task.products)
        return list(task_dependencies - task_products)

    def check_dependencies(self):
        self._logmessage("Checking dependencies")
        deps = self.dependencies()
        if deps:
            ok = True
            for dependency in deps:
                if os.path.isfile(os.path.join(self.workdir, dependency)):
                    self._logmessage( "   Found  {}".format(os.path.join(self.workdir, dependency)))
                else:
                    self._logmessage("*Missing* {}".format(os.path.join(self.workdir, dependency)))
                    ok = False
            if not ok:
                self._logmessage("ERROR: Required file(s) missing")
                sys.exit(1)
        else:
            self._logmessage("Pipeline has no dependencies")
            return True

    def dry_run(self):
        for task in self.tasks:
            print task

    def run(self):
        for task in self.tasks:
            try:
                proc = task.run()
            except DependencyError as err:
                self._logmessage("ERROR: {} missing for {}".format(', '.join(err.missing), task.name))
            if not proc and task.name:
                self._logmessage("Up to date: {}".format(task.name))
            elif task.outfile != subprocess.PIPE:
                self._logmessage("Starting: {}".format(task.name))
                proc.wait()
                self._logmessage("Finished: {}".format(task.name))


def setup_pipeline(pipefile):
    with open(pipefile) as f:
        pipeline = json.loads(f.read())
    parameters = set()
    parameters.update(re.findall(r"{(.*?)}", pipeline["workdir"]))
    for task in pipeline["tasks"]:
        parameters.update(re.findall(r"{(.*?)}", task["command"]))
    with open("config.txt", "w") as outfile:
        outfile.write("{\n")
        outfile.write('\t"pipeline": "{}",\n'.format(pipefile))
        paramlist = list(parameters)
        for parameter in paramlist[:-1]:
            outfile.write('\t"{}": "",\n'.format(parameter))
        outfile.write('\t"{}": ""\n'.format(paramlist[-1]))
        outfile.write("}\n")


if __name__ == "__main__":
    logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)
    if sys.argv[1] == "Setup":
        setup_pipeline(sys.argv[2])
        print "Pipeline configuration created"
        print "Edit config.txt to specify pipeline parameters"
    elif sys.argv[1] == "Test":
        with open("config.txt") as f:
            pipeline = Pipeline(f.read())
        pipeline.dry_run()
    elif sys.argv[1] == "Run":
        with open("config.txt") as f:
            pipeline = Pipeline(f.read())
        pipeline.run()

