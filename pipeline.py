"""A set of wrapper classes to make designing ruffus pipelines easier.

Provides a standard, extensible interface for ruffus pipelines which
allows the user to design a set of standard functions for assigning
file paths and assessing task up-to-date status.

A pipeline is designed by extending PipelineMetadata with the desired
custom file-path templates, and extending the Task
base-class with the desired argument generator, up-to-date
checking function and execution algorithm.

Then a pipeline is put together at the module level (cannot be inside
classes or functions) by creating instances of the Task
subclass with a metadata argument, and then decorating a function
which calls task_instance() with @ruffus.follows([task-names]),
@ruffus.files(task_instace.arg_gen), and
@check_ifuptodate(task_intance.up2date_checker).

Finally, run the pipeline by calling
ruffus.pipeline_run([*args], [**kwargs])

For example, a pipeline module might look like this:
>>> import ruffus
... class Foo_md(PipelineMetadata):
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
... if __name__ == "__main__":
...     run_args, run_kwargs, md_args, md_kwargs = get_params()
...     md = Foo_md(*md_args, **md_kwargs)
...
...     _task1 = Foo_task1(md)
...     @ruffus.follows()
...     @ruffus.files(_task1.args_gen)
...     @functools.wraps(_task1)
...     def task1(*args, **kwargs):
...         _task1(*args, **kwargs)
...
...     _task2 = Foo_task2(md)
...     @ruffus.follows()
...     @ruffus.files(_task2.args_gen)
...     @functools.wraps(_task2)
...     def task2(*args, **kwargs):
...         _task2(*args, **kwargs
...
...     ruffus.run_pipeline(*run_args, **run_kwargs)

TODO: Update the example.  Won't work currently.
"""
import ruffus
import functools

class BaseMetadata(object):
    _path_tmplts = {}
    _reqd_kwargs = []
    
    def __init__(self, **kwargs):
        for key in self._reqd_kwargs:
            assert key in kwargs
            assert kwargs[key] is not None
        for key in kwargs:
            self.__setattr__(key, kwargs[key])
    
    def get_path(self, file_type, **kwargs):
        return self._path_tmplts[file_type] % dict(list(self.__dict__.items()) + \
                                                   list(kwargs.items()))    

class BaseTask(object):
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        
    def args_gen(self, *args, **kwargs):
        raise NotImplementedError("You must define args_gen for each Task subclass.")
    
    def up2date_checker(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define up2date_checker for each Task subclass.")
        # we don't want to make pipeline think that it has an up2date_checker
        # method when it really doesn't.
    
    def __call__(self, inputs_list, outputs_list, *args, **kwargs):
        raise NotImplementedError("You must define __call__ for each Task subclass.")
    
class BasePipeline(object):
    _topology = {}
    
    def __init__(self, pipeline_metadata):
        self.pipeline_metadata = pipeline_metadata
        for task_class in self._topology:
            self._initialize(task_class)
            
    def _function(self, task_class):
        class_string = self.__class__.__name__.lower()
        task_string = task_class.__name__.lower()
        function_string = "{class_string}__{task_string}".format(class_string = \
                                                                 class_string, 
                                                                 task_string = \
                                                                 task_string)
        return function_string
    
    def _dependencies(self, task_class):
        dep_classes = self._topology[task_class]
        dep_functions = []
        for dep_class in dep_classes:
            dep_functions += [self._function(dep_class)]
        return dep_functions
    
    def _initialize(self, task_class):
        if hasattr(task_class, 'up2date_checker'):
            code = """
@ruffus.follows(self._dependencies(task_class))
@ruffus.files(task_object.args_gen)
@ruffus.check_if_uptodate(task_object.up2date_checker)
@functools.wraps(task_class, __doc__)
def {function}(*args, **kwargs):
    return task_object(*args, **kwargs)
"""
        else:
            code = """
@ruffus.follows(self._dependencies(task_class))
@ruffus.files(task_object.args_gen)
@functools.wraps(task_class, '__doc__')
def {function}(*args, **kwargs):
    return task_object(*args, **kwargs)
"""
        code = code.format(function = self._function(task_class))
        print(code)
        exec(code, globals(), dict(self = self,
                                   task_class = task_class,
                                   task_object = task_class(self.pipeline_metadata)))
    
    def __call__(self):
        pass
        
        
        
        
        
        
def test_exec(function_name):
    class FooTask(BaseTask):
        """Docstring"""
        pass
    def d(task_class):
        return ['some', 'depends']
    self = BasePipeline('md')
    self._dependencies = d
    task_class = FooTask
    task_object = FooTask(self.pipeline_metadata)
    code = """
@ruffus.follows(*self._dependencies(task_class))
@ruffus.files(task_object.args_gen)
@ruffus.check_if_uptodate(task_object.up2date_checker)
@functools.wraps(task_class, ['__doc__'])
def {function}(*args, **kwargs):
    return task_object(*args, **kwargs)
"""
    code = code.format(function = function_name)
    exec(code, globals(), dict(ruffus = ruffus,
                               self = self,
                               task_class = task_class,
                               task_object = task_object))
        
        