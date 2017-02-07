import subprocess
import sys
import os
import shlex
import json
import re
import logging
import collections
import multiprocessing

class DependencyError(Exception):
    def __init__(self, missing):
        self.missing = missing

class Task:
    def __init__(self, commandline, name, parent, infile=None, outfile=None, dependencies=None, products=None, workdir=None):
        self.commandline = commandline
        self.name = name
        self.parent = parent
        self.dependencies = dependencies
        self.products = products
        self.workdir = workdir

    def __str__(self):
        return self.commandline

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
                self._logmessage("ERROR: {} missing".format(', '.join(missing)))
                raise DependencyError()

    def _start_process(self, commandline, pipe=None):
        if '|' in commandline:
            left, right = commandline.split('|', 1)
            leftproc = self._start_process(left, pipe=subprocess.PIPE)
            return self._start_process(right, pipe=leftproc.stdout)
        else:
            command, infile, outfile = self._parse_command(commandline)
            stdin = None
            stdout = None
            if pipe == subprocess.PIPE:
                if infile:
                    stdin = open(os.path.join(self.workdir, infile), 'r')
                stdout = subprocess.PIPE
            elif pipe:
                stdin = pipe
                if outfile:
                    stdout = open(os.path.join(self.workdir, outfile), 'wb')
            else:
                if outfile:
                    stdout = open(os.path.join(self.workdir, outfile), 'wb')
                if infile:
                    stdin = open(os.path.join(self.workdir, infile), 'r')
            return subprocess.Popen(shlex.split(command), stdin=stdin, stdout=stdout, cwd=self.workdir)

    def _logmessage(self, message):
        self.parent._logmessage("{} - {}".format(self.name, message))

    def dry_run(self):
        try:
            if self.outdated():
                self._logmessage(self.commandline)
        except OSError:
            self._logmessage(self.commandline)

    def run(self):
        self.check_dependencies()
        if not self.outdated():
            self._logmessage("Up to date")
        else:
            self._logmessage("Starting")
            self._start_process(self.commandline).wait()
            self._logmessage("Finished")


class Pipeline:
    def __init__(self, parameters):
        self.parameters = parameters
        logging.info("Loading pipeline from {}".format(self.parameters["pipeline"]))
        with open(self.parameters["pipeline"]) as f:
            pipeline = json.loads(f.read())
        self.name = pipeline["name"]
        self.steps = pipeline["tasks"]
        if "workdir" in pipeline:
            self.workdir = pipeline["workdir"].format(**self.parameters)
        else:
            self.workdir = ""
        self._initialize()
        self.check_dependencies()

    def _logmessage(self, message):
        if self.workdir:
            logging.info("{} in {}: {}".format(self.name, self.workdir, message))
        else:
            logging.info("{}: {}".format(self.name, message))

    def parameters(self):
        parameters = set()
        for command in self.commands:
            parameters.update(re.findall(r"{(.*?)}", command))
        return list(parameters)

    def _initialize(self):
        self._logmessage("Initializing")
        #Setup tasks
        self.tasks = []
        for step in self.steps:
            if "command" in step:
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
                self.tasks.append(Task(commandline=commandline, name=step["name"], parent=self, dependencies=dependencies, products=products, workdir=self.workdir))
            elif "pipeline" in step:
                params = self.parameters.copy()
                for k, v in self.parameters[step["name"]].items():
                    if v:
                        params[k] = v
                params["pipeline"] = os.path.join(os.path.dirname(self.parameters["pipeline"]), step["pipeline"])
                pipeline = Pipeline(params)
                self.tasks.append(pipeline)
        #Set dependencies and products
        task_dependencies = set()
        task_products = set()
        for task in self.tasks:
            if isinstance(task, Task):
                if task.dependencies:
                    task_dependencies.update(task.dependencies)
                if task.products:
                    task_products.update(task.products)
        self.dependencies = list(task_dependencies - task_products)
        self.products = list(task_products)

    def check_dependencies(self):
        self._logmessage("Checking dependencies")
        deps = self.dependencies
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
        self._logmessage("Dry run started")
        for task in self.tasks:
            task.dry_run()
        self._logmessage("Dry run finished")

    def run(self):
        self._logmessage("Pipeline started")
        for task in self.tasks:
            task.run()
        self._logmessage("Pipeline finished")


def generate_config(pipefile):
    with open(pipefile) as f:
        pipeline = json.loads(f.read())
    paramset = set()
    parameters = {}
    paramcounts = collections.defaultdict(int)
    for k, v in pipeline.items():
        if isinstance(v, basestring): #NOT PYTHON3 compatible. Need to change later
            paramset.update(re.findall(r"{(.*?)}", v))
    for task in pipeline["tasks"]:
        if "command" in task:
            paramset.update(re.findall(r"{(.*?)}", task["command"]))
        elif "pipeline" in task:
            parameters[task["name"]], subpcounts = generate_config(os.path.join(os.path.dirname(pipefile), task["pipeline"]))
            for param in subpcounts.keys():
                paramcounts[param] += 1
    for param in paramset:
        if "prevent_global" in pipeline and param in pipeline["prevent_global"]:
            continue
        paramcounts[param] += 1
    parameters.update({k: "" for k in paramset})
    return parameters, paramcounts


if __name__ == "__main__":
    logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)
    if sys.argv[1] == "Setup":
        if os.path.isfile(sys.argv[3]):
            print "{} already exists. Stopping.".format(sys.argv[3])
        else:
            config, paramcounts = generate_config(sys.argv[2])
            outstr = json.dumps(config, sort_keys=True, indent=4, separators=(',', ': ')).strip()[1:-1] #Remove outer braces
            extra = '    "pipeline": "{}",\n'.format(sys.argv[2])
            for k, v in paramcounts.items():
                if v > 1:
                    extra += '    "{}": "",\n'.format(k)
            with open(sys.argv[3], 'w') as outfile:
                outfile.write("{\n" + extra + outstr + "}\n")
            print "Pipeline configuration created"
            print "Edit {} to specify pipeline parameters".format(sys.argv[3])
    elif sys.argv[1] == "Test":
        with open(sys.argv[2]) as f:
            pipeline = Pipeline(json.loads(f.read()))
        pipeline.dry_run()
    elif sys.argv[1] == "Run":
        with open(sys.argv[2]) as f:
            pipeline = Pipeline(json.loads(f.read()))
        pipeline.run()

