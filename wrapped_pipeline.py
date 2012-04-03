"""

"""
from ruffus_wrapper import *
import ruffus
import functools

class Foo_md(PipelineMetadata):
    _path_tmplts = dict(input = 'test.in', output = 'test.out')

class Foo_task1(Task):
    "Make an input file."
    def args_gen(self):
        yield([], [self.pipeline_metadata.get_path('input')])
        
    def __call__(self, input_path_list, output_path_list):
        with open(output_path_list[0], 'w') as output_file:
            output_file.write("This is\nAn Output\n file.")

class Foo_task2(Task):
    "Copy input file to output file."
    def args_gen(self):
        yield([self.pipeline_metadata.get_path('input')],
              [self.pipeline_metadata.get_path('output')])

    def __call__(self, input_path_list, output_path_list):
        with open(input_path_list[0]) as input_file, \
             open(output_path_list[0], 'w') as output_file:
            for line in input_file:
                output_file.write(line)

if __name__ == "__main__":
    md = Foo_md()

    _task1 = Foo_task1(md)
    @ruffus.follows()
    @ruffus.files(_task1.args_gen)
    @functools.wraps(Foo_task1, assigned = ['__doc__'])
    def task1(*args, **kwargs):
        _task1(*args, **kwargs)

    _task2 = Foo_task2(md)
    @ruffus.follows(task1)
    @ruffus.files(_task2.args_gen)
    @functools.wraps(Foo_task2, assigned = ['__doc__'])
    def task2(*args, **kwargs):
        _task2(*args, **kwargs)

    ruffus.pipeline_run([task2], verbose = 5, multiprocess = 2)
    