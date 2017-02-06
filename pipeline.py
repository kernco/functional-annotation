import subprocess
import sys
import os
import shlex
import json
import re
import logging


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
        self.parameters = json.loads(jsonstr)
        logging.info("Loading pipeline from {}".format(self.parameters["pipeline"]))
        with open(self.parameters["pipeline"]) as f:
            pipeline = json.loads(f.read())
        self.name = pipeline["name"]
        self.steps = pipeline["tasks"]
        #self.commands = [task["command"] for task in pipeline["tasks"]]
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

    def parameters(self):
        parameters = set()
        for command in self.commands:
            parameters.update(re.findall(r"{(.*?)}", command))
        return list(parameters)

    def _initialize(self):
        logging.info("Initializing {}".format(self.name))
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
            self.tasks.append(Task(command, infile=infile, outfile=outfile, dependencies=dependencies, products=products))
            if len(commands) > 1:
                self.tasks[-1].outfile = subprocess.PIPE
                for command, infile, outfile in [self._parse_command(x) for x in commands[1:-1]]:
                    self.tasks.append(Task(command, infile=self.tasks[-1], outfile=subprocess.PIPE, dependencies=dependencies, products=products))
                command, infile, outfile = self._parse_command(commands[-1])
                self.tasks.append(Task(command, infile=self.tasks[-1], outfile=outfile, dependencies=dependencies, products=products))

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
        logging.info("Checking dependencies")
        deps = self.dependencies()
        if deps:
            ok = True
            for dependency in deps:
                line = "{}: ".format(dependency)
                if os.path.isfile(dependency):
                    line += "Found"
                    logging.info(line)
                else:
                    line += "*Not Found*"
                    ok = False
                    logging.error(line)
            if not ok:
                logging.error("ERROR: Required file(s) missing")
                sys.exit(1)
        else:
            logging.info("Pipeline has no dependencies")
            return True

    def dry_run(self):
        for task in self.tasks:
            print task

    def run(self):
        for task in self.tasks:
            logging.info(str(task))
            proc = task.run()
            if task.outfile != subprocess.PIPE:
                proc.wait()


def setup_pipeline(pipefile):
    with open(pipefile) as f:
        pipeline = json.loads(f.read())
    parameters = set()
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
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
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

