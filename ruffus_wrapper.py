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
...     def __call__(self, input_path_list, output_path_list):
...         with open(output_path_list[0], 'w') as output_file:
...             output_file.write("This is\nAn Output\n file.")
...
... class Foo_task2(Task):
...     "Copy input file to output file."
...     def args_gen(self):
...         yield([self.pipeline_metadata.get_path('input')],
...               [self.pipeline_metadata.get_path('output')])
...     def __call__(self, input_path_list, output_path_list):
...         with open(input_path_list[0]) as input_file, open(output_path_list[0], 'w') as output_file:
...             for line in input_file:
...                 output_file.write(line)
...
... class Foo_pipeline(Pipeline):
...     _task_classes = set([Foo_task1, Foo_task2])
...     _topology = dict(Foo_task1 = [], Foo_task2 = [Foo_task1])
...
... md = Foo_md()
... f = Foo_pipeline(md)
... f.run()

"""
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
    
    def up2date_checker(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define up2date_checker for each Task subclass.\n\
If it is not defined the default ruffus up-to-date checker will run.")
        # Consider removing this from the base_class, since I don't think
        # we want to make pipeline think that it has an up2date_checker
        # method when it really doesn't.
    
    def __call__(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define __call__ for each Task subclass.")

class Pipeline(object):
    _task_classes = set()
    _topology = {}
    
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        self._initialized = {}
        self._tasks = []
        for task_class in self._task_classes:
            self._tasks += [task_class(self.pipeline_metadata)]
        for task in self._tasks:
            self.initialize(task)
        
    def get_follows(self, task):
        follows_list = self._topology.get(task.__class__)
        return follows_list
    
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
                @ruffus.follows(*self.get_follows(task))
                @ruffus.files(task.args_gen)
                @ruffus.check_if_uptodate(task.up2date_checker)
                def wrapped_task(*args, **kwargs):
                    task(*args, **kwargs)
            else:
                @ruffus.follows(*self.get_follows(task))
                @ruffus.files(task.args_gen)
                def wrapped_task(*args, **kwargs):
                    task(*args, **kwargs)
            # It looks like I should be hacking a bunch of other
            # attributes in ruffus_nodes. :(
            # TODO:
            method_name = "_%s" % task.__class__.__name__.lower()
            self.__setattr__(method_name, wrapped_task)
            # This whole block here is a hack to change the ruffus-saved name of
            # the task to what I want, which is the Task class-name.
            dummy_node_name = '%s.wrapped_task' % wrapped_task.__module__
            new_node_name = '%s.%s' % (__name__, task.__class__.__name__)
            # is it possible that by using self.__class__ instead of
            # wrapped_task.__module__, that I could define multiple pipelines
            # in one module?  This would be optimal! 
            ruffus_node = ruffus_nodes.lookup_node_from_name(dummy_node_name)
            ruffus_nodes._name_to_node[new_node_name] = ruffus_node
            ruffus_node.semaphore_name = new_node_name
            del ruffus_nodes._name_to_node[dummy_node_name]
            # So will this task (node in ruffus speak) be 'findable' by other tasks with
            # incoming or outgoing dependencies?
            self._initialized[task] = True
            
    def all_initialized(self):
        for task in self.tasks():
            if not self.initialized(task):
                return False
        return True
            
    def run(self, *args, **kwargs):
        ruffus.pipeline_run([self.tasks[-1]])
    
