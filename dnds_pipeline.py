"""

"""
import functools
import optparse
import os
import ruffus
from pipeline import *
from utils import *

DEFAULT_TOP_DIR = '.'

class Metadata(BaseMetadata):
    _path_tmplts = dict()

class MakeTopDir(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        os.makedirs(path = output_path_list[0])
    
    def up2date_checker(self):
        pass
    
class MakeCodemlDir(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        os.makedirs(path = output_path_list[0])
    
    def up2date_checker(self):
        pass
    
class TranslateFn(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        translate_file(fn_path = input_path_list[0], fa_path = output_path_list[0])
    
    def up2date_checker(self):
        pass
    
class AlignFa(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        align_file(fa_path = input_path_list[0], afa_path = output_path_list[0])
    
    def up2date_checker(self):
        pass
    
class BackAlignFn(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        backalign_file(fn_path = input_path_list[0],
                       afa_path = input_path_list[1],
                       afn_path = input_path_list[0])
    
    def up2date_checker(self):
        pass
    
class SubsampleSeqs(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        # TODO: Fill in everything starting here.
        pass
    
    def up2date_checker(self):
        pass
    
class CalcTree(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class PareTree(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class MakeLabels(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class LabelTree(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class Fa2Phy(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class RunCodeml(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class ProcessTestResults(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass
    
class CompileAllResults(BaseTask):
    def args_gen(self):
        pass
        
    def __call__(self, input_path_list, output_path_list):
        pass
    
    def up2date_checker(self):
        pass


def get_user_params():
    usage = "usage: %prog [options] fn-file design-file"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-d", "--dir", dest = "top_dir",
                      default = DEFAULT_TOP_DIR,
                      help = "top level directory for intermediate analysis \
                              files.  DEFAULT: '%s'" % DEFAULT_TOP_DIR)
    parser.add_option("-N", "--name", "--prefix",
                      dest = "prefix", default = "",
                      help = "prefix to apply to all files residing in the \
                              top-level directory")
    parser.add_option("-r", "--reps", dest = "num_reps", type = 'int',
                      default = 1,
                      help = "number of replicate subsamplings.  \
                              DEFAULT: 1")
    parser.add_option("-n", "--seqs", dest = "num_seqs",
                      type = 'int', default = None,
                      help = "number of sequences in the subsample.  \
                              Increased sample size increases both \
                              power and analysis time.  DEFAULT: use \
                              all sequences.")
    parser.add_option("-p", "--threads", dest = "num_processes",
                      type = int, default = 1,
                      help = "maximum concurrent threads to be run.  \
                              DEFAULT: 1")
    parser.add_option("-v", "--verbose", dest = "verbosity",
                      type = 'int', default = 1,
                      help = "how much to print to stderr.  DEFAULT: 1")
    parser.add_option("-q", "--quiet", dest = "verbosity",
                      action = "store_const", const = 0,
                      help = "print nothing to the screen. DEFAULT: False")
    opts, args = parser.parse_args()    
    fn_path = args[0]
    design_path = args[1]
    verbosity = opts.verbosity


if __name__ == "__main__":
    run_args, run_kwargs, metadata_args, metadata_kwargs = get_user_params()
    metadata = Metadata(*metadata_args, **metadata_kwargs)
    # ---DEFINE ALL WRAPPED TASKS HERE---
    # Just copy paste the template and replace [Task] with the
    # subclass of pipeline.Task, [dependencies] with all functions
    # which must be run before [task], and [task] with the name for the
    # function which should call _task(*args, **kwargs).
    #
    #    _task = [Task](metadata)
    #    @ruffus.follows([dependencies])
    #    @ruffus.files(_task.args_gen)
    #    @functools.wraps(_task.__class__, assigned = ['__doc__'])
    #    def [task](*args, **kwargs):
    #        _task(*args, **kwargs)
    # 
    _task = MakeTopDir(metadata)
    @ruffus.follows()
    @ruffus.files(_task.args_gen)
    @functools.wraps(_task.__class__, assigned = ['__doc__'])
    def make_top_dir(*args, **kwargs):
        _task(*args, **kwargs)
    # ---#############################---
    ruffus.pipeline_run()