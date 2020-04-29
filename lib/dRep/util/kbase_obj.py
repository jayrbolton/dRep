import os
import time
import glob
import re
import subprocess
import shutil
import pandas as pd
from Bio import SeqIO, SeqUtils
from functools import reduce
import logging

import drep

from .dprint import dprint
from .config import _globals


class BinnedContigs:
    '''
    DS for BinnedContigs information
    Very mutable, not necessarily in a consistent state

    Instance variables created during init:
    * upa
    * bins_dir
    * name
    * assembly_upa
    * original_bin_name_list

    Other instance variables
    * stats

    '''

    created_instances = [] # loaded from KBase
    saved_instances = [] # saved to KBase
    
      
    def __init__(self, upa, get_bins_dir=None, **kwargs):
       
        self.created_instances.append(self)
        self.upa = upa


        if get_bins_dir == 'download': 
            self._load_bins_dir()
        elif get_bins_dir == 'local': 
            self.bins_dir = kwargs['bins_dir']
        else:
            raise Exception('must specify `get_bins_dir`')

        self._get_obj()
        self.original_bin_name_list = self.get_curr_bin_name_list()


    def _load_bins_dir(self):
        ''''''
        logging.info(f'Downloading files for BinnedContigs with upa: {self.upa}')

        dprint('_globals', run=globals())

        mguObjData = _globals.mgu.binned_contigs_to_file({
            'input_ref': self.upa, 
            'save_to_shock': 0
            }) # returns dict with just bin_file_directory

        self.bins_dir = mguObjData['bin_file_directory']
    
        dprint('os.listdir(self.bins_dir)', run={**locals(), **globals()})


    def _get_obj(self):
        obj = _globals.ws.get_objects2({
            'objects': [{
                'ref': self.upa
                }]
            }) # huge -- includes all the statistics

        self.name = obj['data'][0]['info'][1]
        self.assembly_upa = obj['data'][0]['data']['assembly_ref']

   
    def get_curr_bin_name_list(self):
        bin_name_list = []
        dprint('self.bins_dir', 'os.listdir(self.bins_dir)', run={**globals(), **locals()})
        for bin_name in next(os.walk(self.bins_dir))[2]:
            if not re.search(r'.*\.fasta$', bin_name):
                msg = f'Found non .fasta bin name {bin_name} in dir {self.bins_dir} for BinnedContigs obj {self.name} with UPA {self.upa}'
                log.warning(msg)
            bin_name_list.append(bin_name)
        return bin_name_list


    def get_curr_bin_path_list(self):
        return [os.path.join(self.bins_dir, bin_name) for bin_name in self.get_curr_bin_name_list()]
    

    def is_empty(self):
        return len(self.get_curr_bin_name_list()) == 0


    def save(self):
        ''''''
        if self.is_empty():
            msg ='Not saving empty dereplicated BinnedContigs %s' % self.name 
            logging.warning(msg)
            _globals.warnings.append(msg)
            return

        self.write_reduced_bin_summary()

        objData = _globals.mgu.file_to_binned_contigs({
            'file_directory': self.bins_dir,
            'assembly_ref': self.assembly_upa,
            'binned_contig_name': self.name + '.dRep',
            'workspace_name': _globals.workspace_name
            })

        dprint('objData:', objData, f'for BinnedContigs {self.name}')

        self.saved_instances.append(self)

        return {
            'ref': objData['binned_contig_obj_ref'],
            'description': 'Dereplicated genomes for ' + self.name
        }


    def pool_into(self, binsPooled_dir):
        '''for all bins, modify to unique name and copy into dir_binsPooled''' 
        for bin_name in self.get_curr_bin_name_list():
            bin_name_new = self.transform_bin_name(bin_name)         

            bin_path = os.path.join(self.bins_dir, bin_name)
            bin_path_new = os.path.join(binsPooled_dir, bin_name_new)

            shutil.copyfile(bin_path, bin_path_new) 
        
    
    def transform_bin_name(self, bin_name):
        '''
        Transform bin name to something unique so bins from all BinnedContigs can be pooled
        '''
        return self.upa.replace('/', '-') + '__' + self.name + '__' + bin_name


    def reduce_to_dereplicated(self, bins_derep_dir):
        '''remove bins not in bins_derep_dir'''
        bins_derep_name_list = os.listdir(bins_derep_dir)
        for bin_name in self.get_curr_bin_name_list():
            if self.transform_bin_name(bin_name) not in bins_derep_name_list:
                os.remove(os.path.join(self.bins_dir, bin_name))

       
    def calc_stats(self, use_transformed_name=True):
        '''
        Calc length, N50, GC, in ds like:
        { 
            'bin_stats': {
                'Bin.001.fasta': {
                    'length': 3004,
                    'N50': 37,
                    'GC': 0.22
                },
                'Bin.012.fasta': {
                    'length': 9005,
                    'N50': 50,
                    'GC': 0.54
                },
                ...
            },
            'length': 55000,
            'GC': 0.48
        }   
        Call before dereplication to get stats for all original bins
        '''
        self.stats = {'bin_stats': {}}
        bin_stats = self.stats['bin_stats']

        bin_name_list = self.get_curr_bin_name_list()
        for bin_name in bin_name_list:
            bin_path = os.path.join(self.bins_dir, bin_name) 
            d = {}
            lens = []
            GCs = []
            for seq in SeqIO.parse(bin_path, 'fasta'):
                lens += [len(seq)]
                GCs += [SeqUtils.GC(seq.seq)]
            
            d['GC'] = reduce(lambda x,y: x+y, map(lambda x,y: x*y, lens, GCs), 0) / sum(lens)
            d['length'] = sum(lens)
            d['N50'] = drep.d_filter.calc_n50(bin_path)

            bin_stats[bin_name] = d

        lens = [bin_stats[bin_name]['length'] for bin_name in bin_name_list]
        GCs = [bin_stats[bin_name]['GC'] for bin_name in bin_name_list]

        self.stats['length'] = sum(lens)
        self.stats['GC'] = reduce(lambda x,y: x+y, map(lambda x,y: x*y, lens, GCs), 0) / sum(lens)

        if use_transformed_name:
            for bin_name in bin_name_list:
                self.stats['bin_stats'][self.transform_bin_name(bin_name)] = self.stats['bin_stats'].pop(bin_name)

        return self.stats 


    def write_reduced_bin_summary(self):
        '''
        Assuming that self was loaded and has had its bins in bins_dir reduced,
        build summary of reduced bins from original BinnedContigs' information
        
        *.summary file could be one of below format:

        Bin name                  Abundance  Completeness    Genome size     GC content
        maxbin_output.001.fasta   0.00       97.2%           2690533         52.9
        
        Bin name                  Completeness    Genome size     GC content
        maxbin_output.001.fasta   97.2%           2690533         52.9
        '''
        stats_dict_list = _globals.ws.get_objects2({
            'objects': [{
                'ref': self.upa
                }]
            })['data'][0]['data']['bins']

        columns = ['Bin name', 'Completeness', 'Genome size', 'GC content']
        smmr = pd.DataFrame(columns=columns)

        bins_reduced_name_list = self.get_curr_bin_name_list()

        for stats_dict in stats_dict_list:
            if stats_dict['bid'] in bins_reduced_name_list:
                smmr.loc[len(smmr)] = [stats_dict[key] for key in ['bid', 'cov', 'sum_contig_len', 'gc']]

        smmr['Completeness'] = smmr['Completeness'].apply(lambda x: str(x * 100) + '%')
        smmr['GC content'] = smmr['GC content'].apply(lambda x: x * 100)

        dprint(smmr.to_csv(sep='\t', index=False))

        with open(os.path.join(self.bins_dir, 'MaxBin2.summary'), 'w') as f:
            f.write(smmr.to_csv(sep='\t', index=False))


    @classmethod
    def clear(cls):
        cls.created_instances = []
        cls.saved_instances = []
