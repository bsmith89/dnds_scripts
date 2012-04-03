"""

"""
from ruffus_wrapper import *
import ruffus
import functools

class Metadata(PipelineMetadata):
    _path_tmplts = dict(input = 'test.in', output = 'test.out')

class Task1(Task):
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

    _task1 = Task1(metadata)
    @ruffus.follows()
    @ruffus.files(_task1.args_gen)
    @functools.wraps(Task1, assigned = ['__doc__'])
    def task1(*args, **kwargs):
        _task1(*args, **kwargs)
        
    ruffus.pipeline_run()