# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
import shutil
import subprocess
import pprint
import uuid
import re
import functools

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils

from .util import report
from .util.dprint import dprint
from .util.kbase_obj import BinnedContigs
from .util.config import _globals, reset

#END_HEADER


class dRep:
    '''
    Module Name:
    dRep

    Module Description:
    A KBase module: dRep
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/n1mus/dRep.git"
    GIT_COMMIT_HASH = "05f421114e50a04008d9897cdb20399b23a47b46"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        callback_url = os.environ['SDK_CALLBACK_URL']
        workspace_url = config['workspace-url']
        shared_folder = config['scratch']
        
        self._globals = { # shared by all API-method runs
            'shared_folder': config['scratch'], 
            'ws': Workspace(workspace_url),
            'dfu': DataFileUtil(callback_url),
            'mgu': MetagenomeUtils(callback_url),
            'kbr': KBaseReport(callback_url),
            }

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def dereplicate(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
 
        #BEGIN dereplicate

        dprint('params', run=locals())


        # set up globals ds `_globals` for this API-method run

        reset(_globals) # clear all fields but `debug`
        _globals.update({
            **self._globals,
            'workspace_name': params['workspace_name'],
            'run_dir': os.path.join(self._globals['shared_folder'], str(uuid.uuid4())),
            'warnings': [],
            })

        os.mkdir(_globals.run_dir)




        #
        ##
        ###
        ####
        #####


        def simple_return(msg, report_params_kwargs={}):
            '''
            Return with a simple message
            '''
            report_params = {
                    'message': msg,
                    'workspace_name': _globals.workspace_name,
                    **report_params_kwargs
                    }
            report_info = _globals.kbr.create_extended_report(report_params)
            output = {
                'report_name': report_info['name'],
                'report_ref': report_info['ref'],
            }

            return [output]



        # 
        ##
        ### check unique UPAs
        #### 
        #####

        if len(set(params['genomes_refs'])) < len(params['genomes_refs']):
            msg = 'Removing duplicate input BinnedContigs'
            logging.warning(msg)
            _globals.warnings.append(msg) 
            params['genomes_refs'] = list(dict.fromkeys(params['genomes_refs']))

       



        # 
        ##
        ### load BinnedContigs files (to scratch) and obj data
        ####
        #####
        

        #
        if _globals.debug and params.get('skip_dl'):
            bins_dir_name_l = params['bins_dir_name_l']
            bins_dir_l = [os.path.join(_globals.shared_folder, bins_dir_name) for bins_dir_name in bins_dir_name_l]

            for upa, bins_dir in zip(params['genomes_refs'], bins_dir_l):
                BinnedContigs(upa, get_bins_dir='local', bins_dir=bins_dir)


        # load from KBase 
        else:

            for upa in params['genomes_refs']:
                bc = BinnedContigs(upa, get_bins_dir='download')
                
                if bc.is_empty():
                    # never mention this bc again
                    params['genomes_refs'].remove(upa)
                    BinnedContigs.created_instances.remove(bc)
                

            if len(params['genomes_refs']) == 0:
                msg = 'Sorry, please input at least one non-empty BinnedContigs'
                raise ValueError(msg)



        #
        ##
        ###
        ####
        #####

        for bc in BinnedContigs.created_instances:
            bc.calc_stats()



        #
        ##
        ### pool
        ####
        #####
        
 
        binsPooled_dir = os.path.join(_globals.run_dir, 'binsPooled')
        os.mkdir(binsPooled_dir)

        for bc in BinnedContigs.created_instances:
            bc.pool_into(binsPooled_dir)
        
        dprint("os.listdir(binsPooled_dir)", run={**locals(), **globals()})


        #
        ##
        ### 
        #### params
        #####
       

        # flatten param groups (if any)

        paramGrp_key_l = ['filtering', 'genome_comparison', 'clustering', 'scoring', 'warnings']

        for paramGrp_key in paramGrp_key_l:
            if paramGrp_key in params:
                paramGrp_d = params[paramGrp_key]
                for paramIndiv_key in paramGrp_d:
                    params[paramIndiv_key] = paramGrp_d[paramIndiv_key]
                params.pop(paramGrp_key)


        #

        dRep_params = ['--debug']


        #
        
        if _globals.debug and 'processors' in params:
            dRep_params.extend(['--processors', params.get('processors')])
        else:
            dRep_params.extend(['--processors', '8'])

                
        # extract non-default parameters

        dRep_param_defaults = {
                'length': 50000,
                'completeness': 75,
                'contamination': 25, 
                'ignoreGenomeQuality': "False", 
                'MASH_sketch': 1000,
                'S_algorithm': 'ANImf',
                'n_PRESET': 'normal',
                'P_ani': 0.9,
                'S_ani': 0.99,
                'SkipMash': "False", 
                'SkipSecondary': "False",
                'cov_thresh': 0.1,
                'coverage_method': 'larger',
                'clusterAlg': 'average', 
                'completeness_weight': 1,
                'contamination_weight': 5,
                'strain_heterogeneity_weight': 1, 
                'N50_weight': 0.5,
                'size_weight': 0,
                'warn_dist': 0.25, 
                'warn_sim': 0.98,
                'warn_aln': 0.25,
                'checkM_method': 'lineage_wf',
                }

        params_bool = [key for key in dRep_param_defaults if dRep_param_defaults[key] == 'False']

        for key in dRep_param_defaults:
            if params.get(key) and params[key] != dRep_param_defaults[key]:
                dRep_params.append('--' + key)
                if key not in params_bool:
                    dRep_params.append(params[key])

        dRep_params = [str(param) for param in dRep_params]



        #
        ##
        ### run dRep
        ####
        #####

        if _globals.debug and params.get('skip_dRep'):

            dRep_workDir = os.path.join(_globals.shared_folder, params['dRep_workDir_name'])

            dRep_cmd = ("Skipped running dRep ... "
                        "but dRep_params are: "
                        f"`{' '.join(dRep_params)}` and workDir is `{dRep_workDir}`")

        else:
            dRep_workDir = os.path.join(_globals.run_dir, 'dRep_workDir')

            dRep_cmd = ' '.join(
                [f'dRep dereplicate {dRep_workDir} --genomes {binsPooled_dir}/*.fasta'] + dRep_params)

            logging.info('Running dRep command `%s`' % dRep_cmd)
            completed_proc = subprocess.run(dRep_cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)

            retcode = completed_proc.returncode

            dprint('retcode', run=locals())

            if retcode != 0:

                # check: if Bdb.csv is empty -> nothing passed length/qual filtering
                with open(os.path.join(dRep_workDir, 'data_tables/Bdb.csv')) as f:
                    num_lines = sum(1 for line in f)

                if num_lines == 1:
                    msg = 'Sorry, dRep terminated with return code %d because no bins passed length and CheckM filtering' % retcode

                else:
                    msg = (
                        "Sorry, commmand: `%s` terminated abnormally with return code: %d"
                        % dRep_cmd, retcode
                        )

                raise Exception(msg)





        #
        ##
        ### result BinnedContigs
        ####
        #####

        bins_derep_dir = os.path.join(dRep_workDir, 'dereplicated_genomes')
        objects_created = []

        # for each original BinnedContigs
        for bc in BinnedContigs.created_instances:
            logging.info(f'Dereplicating {bc.name}')
            bc.reduce_to_dereplicated(bins_derep_dir)

            if not bc.is_empty():
                if _globals.debug and params.get('skip_save_bc'):
                    BinnedContigs.saved_instances.append(bc)
                else: 
                    objects_created.append(bc.save())
            



        #
        ##
        ### html
        ####
        #####

        html_dir = os.path.join(_globals.run_dir, 'html_dir')
        shutil.copytree('/kb/module/ui/output', html_dir) # dir of html and accessories

        report.HTMLBuilder(BinnedContigs.saved_instances, dRep_cmd, dRep_workDir, html_dir)



        #
        ##
        ### to shock
        ####
        #####

        if _globals.debug and params.get('skip_save_shock'):
            return


        def dir_to_shock(dir_path, name, description):
            '''
            For regular directories or html directories
            
            name - for regular directories: the name of the flat (zip) file returned to ui
                   for html directories: the name of the html file
            '''
            dfu_fileToShock_ret = _globals.dfu.file_to_shock({
                'file_path': dir_path,
                'make_handle': 0,
                'pack': 'zip',
                })

            return {
                'shock_id': dfu_fileToShock_ret['shock_id'],
                'name': name,
                'description': description
                }


        workDirZip_shockInfo = dir_to_shock(
            dRep_workDir, 'dRep_work_directory', 'Contains intermediate files, logs, results, etc.')

        htmlZip_shockInfo = dir_to_shock(
            html_dir, 'dRep_dereplicate_report.html', 'dRep dereplicate results')



        #
        ##
        ### Report
        ####
        #####

        report_params = {
                'message': 'dRep dereplicate ran successfully',
                'warnings': _globals.warnings,
                'direct_html_link_index': 0,
                'html_links': [htmlZip_shockInfo],
                'file_links': [workDirZip_shockInfo],
                'report_object_name': 'dRep_report_' + str(uuid.uuid4()),
                'workspace_name': params['workspace_name'],
                'objects_created': objects_created
                }

        report_output = _globals.kbr.create_extended_report(report_params)

        dprint('report_output', run=locals())

        output = {
            'report_name': report_output['name'],
            'report_ref': report_output['ref'],
        }


        #END dereplicate

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method dereplicate return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
