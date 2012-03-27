import ruffus
import optparse
import os
import random
import Bio.SeqIO
import Bio.Align.Applications
import Bio.Application
import backalign
import paretree
import treeassign
import Bio.Phylo.PAML.codeml

DEFAULT_TOP_DIR = "."

DESIGN_PATH_TMPL = "master.design"
MSTR_FN_PATH_TMPL = "master.fn"
MSTR_FA_PATH_TMPL = "master.fa"
MSTR_AFA_PATH_TMPL = "master.afa"
MSTR_AFN_PATH_TMPL = "master.afn"
MSTR_TRE_PATH_TMPL = "master.tre"
LABELS_PATH_TMPL = "master.h%(hypo)d.labels"
NAMES_PATH_TMPL = "rep%(rep)02d.names"
AFN_PATH_TMPL = "rep%(rep)02d.afn"
PHY_PATH_TMPL = "rep%(rep)02d.phy"
TRE_PATH_TMPL = "rep%(rep)02d.tre"
NWK_PATH_TMPL = "rep%(rep)02d.h%(hypo)d.nwk"
CTL_PATH_TMPL = "rep%(rep)02d.h%(hypo)d.ctl"
MLC_PATH_TMPL = "rep%(rep)02d.h%(hypo)d.mlc"
RSLTS_PATH_TMPL = "rep%(rep)02d.results.tsv"
MSTR_RSLTS_PATH_TMPL = "master.results.tsv"
CODEML_DIR_TMPL = "rep%(rep)02d-h%(hypo)d"


class PipelineMetadata():
    def __init__(self, **kwargs):
        self.top_dir = kwargs.get('top_dir')
        self.num_reps = kwargs.get('num_reps')
        self.num_seqs = kwargs.get('num_seqs')
                
    def _get_codeml_dir(self, rep, hypo):
        return os.path.join(self.top_dir, CODEML_DIR_TMPL % dict(rep = rep, hypo = hypo))
        
    def _get_design_path(self):
        return os.path.join(self.top_dir, DESIGN_PATH_TMPL)
        
    def _get_master_fn_path(self):
        return os.path.join(self.top_dir, MSTR_FN_PATH_TMPL)
        
    def _get_master_fa_path(self):
        return os.path.join(self.top_dir, MSTR_FA_PATH_TMPL)
        
    def _get_master_afa_path(self):
        return os.path.join(self.top_dir, MSTR_AFA_PATH_TMPL)
          
    def _get_master_afn_path(self):
        return os.path.join(self.top_dir, MSTR_AFN_PATH_TMPL)
        
    def _get_master_tree_path(self):
        return os.path.join(self.top_dir, MSTR_TRE_PATH_TMPL)
        
    def _get_labels_path(self, hypo):
        return os.path.join(self.top_dir, LABELS_PATH_TMPL % dict(hypo = hypo))
    
    def _get_names_path(self, rep):
        return os.path.join(self.top_dir, NAMES_PATH_TMPL % dict(rep = rep))
    
    def _get_afn_path(self, rep):
        return os.path.join(self.top_dir, AFN_PATH_TMPL % dict(rep = rep))
        
    def _get_phy_path(self, rep):
        return os.path.join(self.top_dir, PHY_PATH_TMPL % dict(rep = rep))
    
    def _get_tree_path(self, rep):
        return os.path.join(self.top_dir, TRE_PATH_TMPL % dict(rep = rep))
        
    def _get_nwk_path(self, rep, hypo):
        return os.path.join(self._get_codeml_dir(rep, hypo),
                            NWK_PATH_TMPL % dict(rep = rep, hypo = hypo))
    
    def _get_ctl_path(self, rep, hypo):
        return os.path.join(self._get_codeml_dir(rep, hypo),
                            CTL_PATH_TMPL % dict(rep = rep, hypo = hypo))
    
    def _get_mlc_path(self, rep, hypo):
        return os.path.join(self._get_codeml_dir(rep, hypo),
                            MLC_PATH_TMPL % dict(rep = rep, hypo = hypo))
    
    def _get_results_path(self, rep):
        return os.path.join(self.top_dir, RSLTS_PATH_TMPL % dict(rep = rep))
    
    def _get_master_results_path(self):
        return os.path.join(self.top_dir, MSTR_RSLTS_PATH_TMPL)
    
    def _make_codeml_dir_params(self):
        """Return an iterator of parameters for make_codeml_dir.
        
        A generator to pass as an on-the-fly parameter generator for
        make_codeml_dir.
        
        """
        for rep in range(self.num_reps):
            for hypo in [0, 1]:
                codeml_dir = self._get_codeml_dir(rep, hypo)
                yield [None, codeml_dir]
        
    def _translate_fn_params(self):
        """Return an iterator of parameters for translate_fn.
        
        A generator to pass as an on-the-fly parameter generator for
        translate_fn.
        
        """
        fn_path = self._get_master_fn_path()
        fa_path = self._get_master_fa_path()
        yield [fn_path, fa_path]
        
    def _align_fa_params(self):
        """Return an iterator of parameters for align_fa.
        
        A generator to pass as an on-the-fly parameter generator for
        align_fa.
        
        """
        fa_path = self._get_master_fa_path()
        afa_path = self._get_master_afa_path()
        yield [fa_path, afa_path]
    
    def _backalign_fn_params(self):
        """Return an iterator of parameters for backalign_fn.
        
        A generator to pass as an on-the-fly parameter generator for
        backalign_fn.
        
        """
        fn_path = self._get_master_fn_path()
        afa_path = self._get_master_afa_path()
        afn_path = self._get_master_afn_path()
        yield [[fn_path, afa_path], afn_path]
    
    def _subsample_seqs_params(self):
        """Return an iterator of parameters for backalign_fn.
        
        A generator to pass as an on-the-fly parameter generator for
        subsample_seqs.
        
        """
        design_path = self._get_design_path()
        afn_path = self._get_master_afn_path()
        for rep in range(self.num_reps):
            afn_out_path = self._get_afn_path(rep)
            names_out_path = self._get_names_path(rep)
            yield [[design_path, afn_path], [names_out_path, afn_out_path], self.num_seqs]
    
    def _calc_tree_params(self):
        """Return an iterator of parameters for calc_trees.
        
        A generator to pass as an on-the-fly parameter generator for
        calc_trees.
        
        """
        afa_path = self._get_master_afa_path()
        tree_path = self._get_master_tree_path()
        yield [afa_path, tree_path]
    
    def _pare_tree_params(self):
        """Return an iterator of parameters for pare_tree.
        
        A generator to pass as an on-the-fly parameter generator for
        pare_tree.
        
        """
        master_tree_path = self._get_master_tree_path()
        for rep in range(self.num_reps):
            names_path = self._get_names_path(rep)
            tree_path = self._get_tree_path(rep)
            yield [[master_tree_path, names_path], tree_path]
    
    def _make_labels_params(self):
        """Return an iterator of parameters for make_labels.
        
        A generator to pass as an on-the-fly parameter generator for
        make_labels.
        
        """
        design_path = self._get_design_path()
        h0_labels_path = self._get_labels_path(0)
        h1_labels_path = self._get_labels_path(1)
        yield [design_path, [h0_labels_path, h1_labels_path]]
    
    def _label_tree_params(self):
        """Return an iterator of parameters for label_tree.
        
        A generator to pass as an on-the-fly parameter generator for
        label_tree.
        
        """
        design_path = self._get_design_path()
        for rep in range(self.num_reps):
            tree_path = self._get_tree_path(rep)
            for hypo in [0, 1]:
                labels_path = self._get_labels_path(hypo)
                nwk_path = self._get_nwk_path(rep, hypo)
                yield [[tree_path, design_path, labels_path], nwk_path]
        
    def _fa2phy_params(self):
        """Return an iterator of parameters for fa2phy.
        
        A generator to pass as an on-the-fly parameter generator for
        fa2phy.
        
        """
        for rep in range(self.num_reps):
            afn_path = self._get_afn_path(rep)
            phy_path = self._get_phy_path(rep)
            yield [afn_path, phy_path]

    def _run_codeml_params(self):
        """Return an iterator of parameters for run_codeml.
        
        A generator to pass as an on-the-fly parameter generator for
        run_codeml.
        
        """
        for rep in range(self.num_reps):
            phy_path = self._get_phy_path(rep)
            for hypo in [0, 1]:
                working_dir = self._get_codeml_dir(rep, hypo)
                nwk_path = self._get_nwk_path(rep, hypo)
                mlc_path = self._get_mlc_path(rep, hypo)
                yield [[working_dir, phy_path, nwk_path], mlc_path]
    
    def _process_results_params(self):
        """Return an iterator of parameters for process_results.
        
        A generator to pass as an on-the-fly parameter generator for
        process_results.
        
        """
        for rep in range(self.num_reps):
            results_path = self._get_results_path(rep)
            mlc_paths = []
            for hypo in [0, 1]:
                mlc_paths += [self._get_mlc_path(rep, hypo)]
            yield [mlc_paths, results_path]
    
    def _combine_results_params(self):
        """Return an iterator of parameters for combine_results.
        
        A generator to pass as an on-the-fly parameter generator for
        combine_results.
        
        """
        master_results_path = self._get_master_results_path()
        results_paths = []
        for rep in range(self.num_reps):
            results_paths += [self._get_results_path(rep)]
        yield [results_paths, master_results_path]

if __name__ == "__main__":
    # TODO:
    # so this isn't the best way to do this.  Instead I should define all of
    # the functions externally from this if statement, and then just decorate
    # them here.
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage = usage)
    parser.add_option("-d", "--dir", dest = "top_dir",
                      default = DEFAULT_TOP_DIR,
                      help = "top level directory for intermediate analysis \
                              files.  DEFAULT: '%s'" % DEFAULT_TOP_DIR)
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
    opts = parser.parse_args()[0]    
    
    
    pipeline_metadata = PipelineMetadata(top_dir = opts.top_dir,
                                         num_reps = opts.num_reps,
                                         num_seqs = opts.num_seqs)
    
    @ruffus.files(pipeline_metadata._make_codeml_dir_params)         
    def _make_codeml_dir(null_input, codeml_dir):
        try:
            os.makedirs(codeml_dir)
        except OSError:
            pass
    
    @ruffus.files(pipeline_metadata._translate_fn_params)
    def _translate_fn(fn_path, fa_out_path):
        out_recs = []
        for rec in Bio.SeqIO.parse(fn_path, 'fasta'):
            rec.seq = rec.seq.ungap(gap = '-').translate()
            out_recs += [rec]
        Bio.SeqIO.write(out_recs, fa_out_path, 'fasta')
        
    @ruffus.follows(_translate_fn)
    @ruffus.files(pipeline_metadata._align_fa_params)
    def _align_fa(fa_path, afa_out_path):
        mafft_call = Bio.Align.Applications.MafftCommandline(input = fa_path)
        mafft_call.maxiterate = 1000
        mafft_call.retree = 2
        stdout = mafft_call()[0]
        with open(afa_out_path, "w") as afa_out_file:
            afa_out_file.write(stdout)
            
    @ruffus.follows(_translate_fn, _align_fa)
    @ruffus.files(pipeline_metadata._backalign_fn_params)
    def _backalign_fn(fn_afa_paths_list, afn_out_path):
        fn_path = fn_afa_paths_list[0]
        afa_path = fn_afa_paths_list[1]
        backalign.files(fn_path, afa_path, afn_out_path)
    
    @ruffus.follows(_backalign_fn)
    @ruffus.files(pipeline_metadata._subsample_seqs_params)
    def _subsample_seqs(design_afn_paths_list, names_afn_out_paths_list, num_seqs):
        design_path = design_afn_paths_list[0]
        afn_path = design_afn_paths_list[1]
        names_out_path = names_afn_out_paths_list[0]
        afn_out_path = names_afn_out_paths_list[1]
        by_treatment = {}
        with open(design_path) as design_file:
            for line in design_file:
                name, treatment = line.strip().split()
                if treatment not in by_treatment:
                    by_treatment[treatment] = []
                by_treatment[treatment] += [name]
        out_seq_names = []
        for treatment in by_treatment:
            treatment_sample = random.sample(by_treatment[treatment], num_seqs)
            out_seq_names.extend(treatment_sample)
        out_recs = []
        with open(names_out_path, 'w') as names_out_file:
            for rec in Bio.SeqIO.parse(afn_path, 'fasta'):
                if rec.name in out_seq_names:
                    names_out_file.write("%s\n" % rec.name)
                    out_recs += [rec]
        assert len(out_recs) == len(out_seq_names)
        Bio.SeqIO.write(out_recs, afn_out_path, 'fasta')
    
    @ruffus.follows(_align_fa)
    @ruffus.files(pipeline_metadata._calc_tree_params)
    def _calc_tree(afa_path, tree_out_path):
        class FastTreeCommandline(Bio.Application.AbstractCommandline):
            def __init__(self, cmd="FastTreeMP", **kwargs):
                self.parameters = []
                self.parameters += \
                    [Bio.Application._Switch(['-nosupport', 'nosupport'],
                             "don't include support values in the output tree")]
                self.parameters += \
                    [Bio.Application._Argument(['', 'input'],
                               'input file')]
                # set a list of parameters which are objects derived from the base class
                # _AbstractParameter
                Bio.Application.AbstractCommandline.__init__(self, cmd, **kwargs)
        fasttree_call = FastTreeCommandline(input = afa_path)
        fasttree_call.nosupport = True
        stdout = fasttree_call()[0]
        open(tree_out_path, 'w').write(stdout)
        
    @ruffus.follows(_calc_tree, _subsample_seqs)
    @ruffus.files(pipeline_metadata._pare_tree_params)
    def _pare_tree(master_tree_names_paths_list, tree_out_path):
        master_tree_path = master_tree_names_paths_list[0]
        names_path = master_tree_names_paths_list[1]
        names = []
        with open(names_path) as names_file:
            for line in names_file:
                names += [line.strip()]
        paretree.files(master_tree_path, names, tree_out_path)
        
    @ruffus.files(pipeline_metadata._make_labels_params)
    def _make_labels(design_path, labels_paths_list):
        h0_labels_out_path = labels_paths_list[0]
        h1_labels_out_path = labels_paths_list[1]
        design = {}
        with open(design_path) as design_file:
            for line in design_file:
                name, treatment = line.strip().split()
                if treatment not in design:
                    design[treatment] = []
                design[treatment] += [name]
        with open(h0_labels_out_path, 'w') as h0_labels_file, \
             open(h1_labels_out_path, 'w') as h1_labels_file:
            i = 1
            for treatment in design:
                h0_labels_file.write("%s\t#1\n" % (treatment))
                h1_labels_file.write("%s\t#%d\n" % (treatment, i))
                i += 1
                
    @ruffus.follows(_make_labels, _pare_tree, _make_codeml_dir)
    @ruffus.files(pipeline_metadata._label_tree_params)
    def _label_tree(tree_design_labels_paths_list, nwk_out_path):
        tree_path = tree_design_labels_paths_list[0]
        design_path = tree_design_labels_paths_list[1]
        labels_path = tree_design_labels_paths_list[2]
        treeassign.files(tree_path, design_path, labels_path, nwk_out_path)
        
    @ruffus.follows(_subsample_seqs)
    @ruffus.files(pipeline_metadata._fa2phy_params)
    def _fa2phy(afn_path, phy_out_path):
        Bio.SeqIO.convert(afn_path, 'fasta', phy_out_path, 'phylip')
        with open(phy_out_path) as phy_file:
            phy_lines = phy_file.readlines()
        first = True
        with open(phy_out_path, 'w') as phy_file:
            for line in phy_lines:
                if first:
                    phy_file.write("%s I\n" % line.strip())
                    first = False
                else:
                    phy_file.write(line)
        
    @ruffus.follows(_fa2phy, _label_tree, _make_codeml_dir)
    @ruffus.files(pipeline_metadata._run_codeml_params)
    def _run_codeml(dir_phy_nwk_paths_list, mlc_path):
        working_dir = dir_phy_nwk_paths_list[0]
        phy_path = dir_phy_nwk_paths_list[1]
        nwk_path = dir_phy_nwk_paths_list[2]
        cml = Bio.Phylo.PAML.codeml.Codeml(alignment = phy_path,
                                           tree = nwk_path,
                                           working_dir = working_dir,
                                           out_file = mlc_path)
        cml.set_options(noisy = 3,
                        verbose = 1,
                        runmode = 0,
                        seqtype = 1,
                        CodonFreq = 2,
                        ndata = 1,
                        clock = 0,
                        aaDist = 0,
                        model = 2,
                        NSsites = [0],
                        icode = 0,
                        Mgene = 0,
                        fix_kappa = 0,
                        kappa = 3.53330,
                        fix_omega = 0,
                        omega = 1,
                        getSE = 0,
                        RateAncestor = 0,
                        Small_Diff = 0.0000005,
                        cleandata = 0,
                        fix_blength = 1,
                        method = 0)
        cml.run(parse = False)
        
    @ruffus.follows(_run_codeml)
    @ruffus.files(pipeline_metadata._process_results_params)
    def _process_results(rep_mlc_paths, rep_results_out_path):
        with open(rep_results_out_path, 'w') as out_file:
            pass
        
    @ruffus.follows(_process_results)
    @ruffus.files(pipeline_metadata._combine_results_params)
    def _report_results(results_paths, master_results_path):
        with open(master_results_path, 'w') as out_file:
            pass
        
    ruffus.pipeline_run([_report_results], multiprocess = opts.num_processes)
    