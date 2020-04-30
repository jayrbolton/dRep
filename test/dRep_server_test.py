# -*- coding: utf-8 -*-
import os
import shutil
import time
import unittest
from configparser import ConfigParser
import functools
import sys
import subprocess
import re
import tarfile
import logging
import json

from installed_clients.WorkspaceClient import Workspace

from dRep.dRepServer import MethodContext
from dRep.authclient import KBaseAuth as _KBaseAuth

from dRep.dRepImpl import dRep
from dRep.util.kbase_obj import BinnedContigs
from dRep.util.dprint import dprint, where_am_i
from dRep.util.config import _globals


SURF_B_2binners = ['34837/23/1', '34837/3/1'] # maxbin, metabat
SURF_B_2binners_CheckM = ['34837/16/1', '34837/2/1', ] # maxbin, metabat
SURF_B_2binners_CheckM_dRep = ['34837/17/13', '34837/18/13'] # maxbin, metabat
capybaraGut_MaxBin2 = ['37096/11/1']
capybaraGut_MetaBAT2 = ['37096/9/1']
capybaraGut_2binners = capybaraGut_MetaBAT2 + capybaraGut_MaxBin2
small_arctic_metabat = ['34837/46/1']


file_combos = {
    'SURF_B_2binners': {
        'genomes_refs': ['34837/23/1', '34837/3/1'], # maxbin, metabat
        'bins_dir_name_l': ['SURF-B.MEGAHIT.maxbin', 'SURF-B.MEGAHIT.metabat']
        },
    'SURF_B_2binners_CheckM': {
        'genomes_refs': ['34837/16/1', '34837/2/1'], # maxbin, metabat
        'bins_dir_name_l': ['SURF-B.MEGAHIT.maxbin.CheckM', 'SURF-B.MEGAHIT.metabat.CheckM'],
        'dRep_workDir_name': 'dRep_workDir_SURF-B.MEGAHIT.2binners.CheckM_taxwf'
        },
    'SURF_B_2binners_CheckM_dRep': {
        'genomes_refs': ['34837/17/13', '34837/18/13'] # maxbin, metabat
        },
    'capybaraGut_MaxBin2': {
        'genomes_refs': ['37096/11/1'],
        'bins_dir_name_l': ['capybaraGut.MaxBin2']
        },
    'capybaraGut_MetaBAT2': {
        'genomes_refs': ['37096/9/1'],
        'bins_dir_name_l': ['capybaraGut.MetaBAT2']
        },
    'capybaraGut_2binners': {
        'genomes_refs': ['37096/11/1', '37096/9/1'],
        'bins_dir_name_l': ['capybaraGut.MaxBin2', 'capybaraGut.MetaBAT2']
        },
    }

param_sets = {

    'numbers': {
        'completeness_weight': 0.5,
        'contamination_weight': 3.5,
        'strain_heterogeneity_weight': 0.6,
        'N50_weight': 0.111,
        'size_weight': -0.3,

        'MASH_sketch': 800,
        'P_ani': 0.8,
        'S_ani': 0.95,
        'cov_thresh': 0.01,

        'percent': 70,

        'warn_dist': 0.255,
        'warn_sim': 1.0,
        'warn_aln': 0.22
        },

    'passive_options': {
        'S_algorithm': 'ANIn',
        'n_PRESET': 'tight',
        'coverage_method': 'total',
        'clusterAlg': 'ward',
        'tax_method': 'max'
        },

    'filtering': {
        'length': 2500000,
        'completeness': 95,
        'contamination': 5
        },


    'ignoreGenomeQuality': {
        'ignoreGenomeQuality': 'True',
        'completeness': 99,
        },
 
    'go_ANI': { # FAIL
        'S_algorithm': 'goANI'
        },


   'ANImf_tight': { # normal may not apply to (default) ANImf
        'S_algorithm': 'ANImf',
        'n_PRESET': 'tight'
        },

   'ANIn_tight': {
       'S_algorithm': 'ANIn',
       'n_PRESET': 'tight'
       },

    'gANI_total': { # total does not apply to gANI
        'S_algorithm': 'gANI',
        'coverage_method': 'total',
        },

    'SkipMASH':  {
        'SkipMASH': 'True',
        'MASH_sketch': 900
        },


    'SkipSecondary': {
        'SkipSecondary': 'True'
        },

    'SkipClustering': { # ?
        'SkipMASH': 'True',
        'SkipSecondary': 'True',
        },

    'SkipAll': { # ?
        'ignoreGenomeQuality': 'True',
        'SkipMASH': 'True',
        'SkipSecondary': 'True',
        }, 
}

#param_sets = {'filtering': param_sets['filtering']}
#param_sets = {key: param_sets[key] for key in list(param_sets.keys())[:-2]}


params_local = {
#---------------- machine specific----------------------    
    'processors': 8,
    'checkM_method': 'taxonomy_wf', # default uses 40GB memory, this uses 16GB (?)
#----------------- skip --------------------------------    
    'skip_dl' : True,
    #'skip_dRep': True,
    'skip_save_bc': True,
    'skip_save_shock': True,
#----------------- test data for skips -----------------    
    **file_combos['SURF_B_2binners_CheckM']
}


class dRepTest(unittest.TestCase):


    # TODO if this throws, kill test suite
    def _test(self):
        ret = self.serviceImpl.dereplicate(self.ctx, {
            'workspace_name': self.wsName,
            'genomes_refs': SURF_B_2binners_CheckM,
            **params_local
            }
        )

    @classmethod
    @where_am_i
    def setUpClass(cls):
        '''
        Run once, after all dRepTest objs have been instantiated
        '''
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser() # does not handle inline comments! or endofline spaces
        config.read(config_file)
        for nameval in config.items('dRep'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'dRep',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = dRep(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.testData_dir = '/kb/module/test/data'
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

        # decompress tarballs
        tarball_l = [os.path.join(cls.testData_dir, f) for f in os.listdir(cls.testData_dir) if re.match(r'^.+\.tar\.gz$', f)]
        dprint('tarball_l', run=locals())
        for tarball in tarball_l:
            tar = tarfile.open(tarball)
            tar.extractall(path=cls.testData_dir)
            tar.close()
        dprint('ls /kb/module/test/data', run='cli')
        cls.listTests()


    @classmethod
    @where_am_i
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')
        tag = ' ' + ('!!!!!!!!!!!!!!!!!!!!!!!!!!' * 40) + ' '
        dprint(tag + 'DO NOT FORGET TO GRAB HTML(S)' + tag)
        cls.listTests()
       
    @classmethod
    def listTests(cls):
        dprint("[name for name, func in cls.__dict__.items() if name.startswith('test') and callable(func)]", run=locals())

    def setUp(self):
        dprint(f"cp -r {self.testData_dir}/* {self.scratch}", run='cli')
        self.listTests()


    def tearDown(self):
        # clear test bins and work dirs
        dprint(f"rm -rf {os.path.join(self.scratch, 'SURF-B*')}", run='cli')
        dprint(f"rm -rf {os.path.join(self.scratch, 'capybara*')}", run='cli')

        #
        BinnedContigs.clear()
        self.listTests()

###################### full network test ###########################################################

    def test_network(self):
        '''
        Often the dl and ul are skipped
        This exercises that code
        '''
        self.serviceImpl.dereplicate(self.ctx, {
            'workspace_name': self.wsName,
            'genomes_refs': SURF_B_2binners_CheckM,
            'processors': 8,
            'ignoreGenomeQuality': 'True',
            }
        )


####################### faulty input tests #########################################################

    def test_dup_BinnedContigs(self):
        self.serviceImpl.dereplicate(self.ctx, {
            'workspace_name': self.wsName,
            **file_combos['SURF_B_2binners_CheckM'],
            'skip_dl' : True,
            'skip_dRep': True,
            'skip_save_bc': True,
            'skip_save_shock': True,
            'genomes_refs': SURF_B_2binners_CheckM * 2 # order matters when `skip_dl`
            }
        )

        self.assertTrue('Removing duplicate input BinnedContigs' in _globals.warnings) # TODO warning/error library to reduce hardcoding?

   
    def test_nothing_passes_filtering(self):
        with self.assertRaises(Exception) as cm:
            self.serviceImple.dereplicate(self.ctx, {
                'workspace_name': self.wsName,
                'genomes_refs': small_arctic_metabat,
                'checkM_method': 'taxonomy_wf',
                }
             )
 
            self.assertTrue(
                'no bins passed length and CheckM filtering' in str(cm.exception))


############################ param combo tests #####################################################

def _gen_test_param_set(params_dRep):
    def test_param_set(self):
        ret = self.serviceImpl.dereplicate(
            self.ctx, 
            {
                'workspace_name': self.wsName,
                **params_dRep,
                **params_local,
            })
    return test_param_set

for (param_set_name, param_set), count in zip(param_sets.items(), range(len(param_sets))):
    test_name = 'test_param_set_' + str(count) + '_' + param_set_name
    #setattr(dRepTest, test_name, _gen_test_param_set(param_set))



########################### decorate test* funcs ###################################################
"""
for key, value in dRepTest.__dict__.items():
    if type(key) == str and key.startswith('test') and callable(value):
        dprint('key', 'value', run=locals())
        setattr(dRepTest, key, where_am_i(value))
"""


############################ select what to run ####################################################
run_tests = ['test_network']

for key, value in dRepTest.__dict__.copy().items():
    if type(key) == str and key.startswith('test') and callable(value):
        dprint('key', 'value', run=locals())
        if key not in run_tests:
            delattr(dRepTest, key)







