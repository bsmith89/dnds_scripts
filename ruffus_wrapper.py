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
... class Foo_task1(Task):
...     "Make an input file."
...     def args_gen(self):
...         yield([], [self.pipeline_metadata.get_path('input')])
...
...     def __call__(self, input_path_list, output_path_list):
...         with open(output_path_list[0], 'w') as output_file:
...             output_file.write("This is\nAn Output\n file.")
...
... class Foo_task2(Task):
...     "Copy input file to output file."
...     def args_gen(self):
...         yield([self.pipeline_metadata.get_path('input')],
...               [self.pipeline_metadata.get_path('output')])
...
...     def __call__(self, input_path_list, output_path_list):
...         with open(input_path_list[0]) as input_file, \
...              open(output_path_list[0], 'w') as output_file:
...             for line in input_file:
...                 output_file.write(line)
...
... class Foo_pipeline(Pipeline):
...     _topology = dict(Foo_task1 = [], Foo_task2 = [Foo_task1])
...
... md = Foo_md()
... f = Foo_pipeline(md)
... f.run()

"""
import functools
import re
import ruffus
from ruffus.graph import node as ruffus_nodes


class PipelineMetadata(object):
    _path_tmplts = {}
    _reqd_kwargs = []
    
    def __init__(self, **kwargs):
        for arg in self._reqd_kwargs:
            assert arg in kwargs
        for key in kwargs:
            self.__setattr__(key, kwargs[key])
    
    def get_path(self, file_type, **kwargs):
        return self._path_tmplts[file_type] % dict(list(self.__dict__.items()) + \
                                                   list(kwargs.items()))    

class Task(object):
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        
    def args_gen(self, *args, **kwargs):
        raise NotImplementedError("You must define args_gen for each Task subclass.")
    
#    def up2date_checker(self, inputs_list, outputs_list, *args, **kwargs):
#        raise NotImplementedError("You must define up2date_checker for each Task subclass.")
#        # we don't want to make pipeline think that it has an up2date_checker
#        # method when it really doesn't.
    
    def __call__(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define __call__ for each Task subclass.")

class Pipeline(object):
    _topology = {}
    
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        # task is always the task-class
        # use self.function() to get the name of the
        # wrapped ruffus function.
        for task in self._tasks():
            self._initialize(task(self.pipeline_metadata))
                
    def _function(self, task):
        """Return the name of the global function which executes the task."""
        return "%s_%s" % (self.__class__.__name__.lower(), task.__name__.lower())
        
    def _get_follows(self, task):
        """Return the names of all functions which task follows."""
        task_list = self._topology.get(task.__class__)
        function_list = []
        for task in task_list:
            function_list += [self._function(task)]
        return function_list
    
    def _tasks(self):
        for key in self._topology:
            yield key
        
    def _initialize(self, task_object):
        if hasattr(task_object, 'up2date_checker'):
            @ruffus.follows(*self._get_follows(task_object))
            @ruffus.files(task_object.args_gen)
            @ruffus.check_if_uptodate(task_object.up2date_checker)
            @functools.wraps(task_object) # what does it need to update?
            def wrapped_task(*args, **kwargs):
                task_object(*args, **kwargs)
        else:
            @ruffus.follows(*self._get_follows(task_object))
            @ruffus.files(task_object.args_gen)
            @functools.wraps(task_object)
            def wrapped_task(*args, **kwargs):
                task_object(*args, **kwargs)
                    
    def __call__(self, *args, **kwargs):
        ruffus.pipeline_run(*args, **kwargs)
    
def reset_ruffus():
    ruffus_nodes._all_nodes = []
    ruffus_nodes._name_to_node = {}
    ruffus_nodes._global_node_index = 0
        
if __name__ == "__main__":
    class Foo_md(PipelineMetadata):
        _path_tmplts = dict(input = 'test.in', output = 'test.out')
    
    class Foo_task1(Task):
        def args_gen(self):
            yield([], [self.pipeline_metadata.get_path('input')])
        
        def __call__(self, inputs_list, outputs_list):
            with open(outputs_list[0], 'w') as outfile:
                outfile.write("This is\nAn output\n f.")
    
    class Foo_task2(Task):
        def args_gen(self):
            yield ([self.pipeline_metadata.get_path('input')],
                   [self.pipeline_metadata.get_path('output')])
        
        def __call__(self, inputs_list, outputs_list):
            with open(inputs_list[0]) as infile, open(inputs_list[0], 'w') as outfile:
                for line in infile:
                    outfile.write(line)
                    
    class Foo_pipeline(Pipeline):
        _topology = {Foo_task1:[], Foo_task2:[Foo_task1]}
    
    f = Foo_pipeline(Foo_md())
    print(ruffus_nodes._all_nodes)
    f()
    reset_ruffus()


















