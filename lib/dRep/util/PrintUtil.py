import functools
import pprint, json
import subprocess

MAX_LINES = 70
print = functools.partial(print, flush=True)

def dprint(*args, **kwargs):
    print('##############################################################')
    for arg in args:
        if isinstance(arg, str):
            print(arg, end=' ', **kwargs)
        else:
            try:
                arg_json = json.dumps(arg,indent=3)
                if arg_json.count('\n') > MAX_LINES:
                    arg_json = '\n'.join(arg_json.split('\n')[0:MAX_LINES] + ['...'])
                print()
                print(arg_json, **kwargs)
            except:
                print(arg, end=' ', **kwargs)
    print()
    print('--------------------------------------------------------------')


def dprint_run(*args, **kwargs):
    for arg in args:
        dprint(arg, **kwargs)
        dprint(subprocess.run(arg, shell=True, stdout=subprocess.PIPE).decode('utf-8'), **kwargs)





