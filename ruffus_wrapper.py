"""A set of wrapper classes for ruffus pipelines.

Provides a standard, extensible interface for ruffus pipelines
which allows the user to design a set of standard functions
for assigning file paths and assessing task up-to-date status.

A pipeline is designed by extending Pipeline with a set of tasks,
while also extending PipelineMetadata with the associated custom
file-path templates, argument generators, and up-to-date status
checkers.

Use like so:
>>> class Foo_md(PipelineMetadata):
...     _path_tmplts = dict(input = 'test.in', output = 'test.out')
...     
... class Foo_task(Task):
...     "Copy input file to output file."
...     def args_gen(self):
...         yield([self.pipeline_metadata.get_path('input')],
...               [self.pipeline_metadata.get_path('output')])
...
...     def __call__(self, input_path_list, output_path_list):
...         with open(input_path_list[0]) as input_file, open(output_path_list[0], 'w') as output_file:
...             for line in input_file:
...                 output_file.write(line)
...
... class Foo_pipeline(Pipeline):
...     _tasks = [Foo_task]
...     _topology = dict(Foo_task = [])
...
... md = Foo_md()
... f = Foo_pipeline(Foo_md)
... f.run()

"""
import ruffus

class PipelineMetadata(object):
    _path_tmplts = {}
    _reqd_kwargs = []
    
    def __init__(self, **kwargs):
        for arg in self._reqd_kwargs:
            assert arg in kwargs
        for key in kwargs:
            self.__setattr__(key, kwargs[key])
    
    def get_path(self, file_type, **kwargs):
        return self._path_tmplts[file_type] % dict(list(self.__dict__.items()),
                                                   list(kwargs.items()))    

class Task(object):
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        
    def args_gen(self, *args, **kwargs):
        raise NotImplementedError("You must define args_gen for each Task subclass.")
    
    def up2date_checker(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define up2date_checker for each Task subclass.\n\
If it is not defined the default ruffus up-to-date checker will run.")
    
    def __call__(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define __call__ for each Task subclass.")

class Pipeline(object):
    _tasks = []
    _topology = {}
    
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        self._initialized = {}
        for i in range(len(self._tasks)):
            # replace a list of Task subclasses with a list of Task objects
            self._tasks[i] = self._tasks[i](self.pipeline_metadata) 
        for task in self._tasks:
            self.initialize(task)
        
    def get_follows(self, task):
        return self._topology[task]
    
    def tasks(self):
        for task in self._tasks:
            yield task
            
    def initialized(self, task):
        if self._initialized.get(task) is True:
            return True
        else:
            return False
        
    def initialize(self, task):
        if not self.initialized(task):
            if hasattr(task, 'up2date_checker'):
                @ruffus.follows(self.get_follows(task))
                @ruffus.files(task.args_gen)
                @ruffus.check_if_uptodate(task.up2date_checker)
                def wrapped_task(*args, **kwargs):
                    task(*args, **kwargs)
            else:
                @ruffus.follows(self.get_follows(task))
                @ruffus.files(task.args_gen)
                def wrapped_task(*args, **kwargs):
                    task(*args, **kwargs)
            self._initialized[task] = True
            
    def all_initialized(self):
        for task in self.tasks():
            if not self.initialized(task):
                return False
        return True
            
    def run(self, *args, **kwargs):
        pass
    
