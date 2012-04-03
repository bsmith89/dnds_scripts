"""

"""
from ruffus_wrapper import *
import ruffus
import functools

class Metadata(PipelineMetadata):
    _path_tmplts = dict(input = 'test.in', output = 'test.out')

class task1(Task):
    "Make an input file."
    def args_gen(self):
        yield([], [self.pipeline_metadata.get_path('input')])
        
    def __call__(self, input_path_list, output_path_list):
        with open(output_path_list[0], 'w') as output_file:
            output_file.write("This is\nAn Output\n file.")


def get_params():
    pass


if __name__ == "__main__":
    run_args, run_kwargs, metadata_args, metadata_kwargs = get_params()
    metadata = Metadata(*metadata_args, **metadata_kwargs)

    _task2 = Foo_task2(md)
    @ruffus.follows(task1)
    @ruffus.files(_task2.args_gen)
    @functools.wraps(Foo_task2, assigned = ['__doc__'])
    def task2(*args, **kwargs):
        _task2(*args, **kwargs)

    ruffus.pipeline_run([task2], verbose = 5, multiprocess = 2)