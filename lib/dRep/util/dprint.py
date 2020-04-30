import functools
import pprint, json
import subprocess
import sys
import os
import time
import logging
import inspect
import time as _time

from .config import _globals


subproc_run = functools.partial(
        subprocess.run, stdout=sys.stdout, stderr=sys.stderr, shell=True, executable='/bin/bash')

TAG_WIDTH = 100
MAX_LINES = 70


def dprint(*args, run=False, where=True, time=False, max_lines=MAX_LINES, subproc_run_kwargs={}, print_kwargs={}):

    if not _globals.debug:
        return

    print = functools.partial(__builtins__['print'], **print_kwargs)

    def print_format(arg):
        if isinstance(arg, (dict, list)):
            arg_json = json.dumps(arg, indent=3, default=str)
            if arg_json.count('\n') > max_lines:
                arg_json = '\n'.join(arg_json.split('\n')[0:max_lines] + ['...'])
            print(arg_json)
        else:
            print(arg)

    print('#' * TAG_WIDTH)

    if where:
        last_frame = inspect.stack()
        print("(file `%s`)\n(func `%s`) " % (os.path.basename(last_frame[1][1]), last_frame[1][3]))
    
    for arg in args:
        if time:
            t0 = _time.time()
        if run:
            print('>> ' + arg)
            if run in ['cli', 'shell']:
                completed_proc = subproc_run(arg, **subproc_run_kwargs)
                retcode = completed_proc.returncode
            elif isinstance(run, dict):
                print_format(eval(arg, run))
            else:
                assert False
        else:
            print_format(arg)
        if time:
            t = _time.time() - t0
            print('[%fs]' % t)
    
    print('-' * TAG_WIDTH)
    
    # return last retcode
    if run and run in ['cli', 'shell']:
        return retcode


def where_am_i(f):
    '''Decorator'''
    def f_new(*args, **kwargs):
        dprint("where am i?\n(func `%s`)" % (f.__qualname__), where=False)
        f(*args, **kwargs)
    return f_new
