from dotmap import DotMap

DEBUG = True # toggle for debugging prints/conditionals

_globals = DotMap(debug=DEBUG) # import for global variables ds

def reset(dm: DotMap):
    dm.clear()
    dm.update({'debug': DEBUG})
