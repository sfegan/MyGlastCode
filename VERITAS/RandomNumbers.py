# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_RandomNumbers', [dirname(__file__)])
        except ImportError:
            import _RandomNumbers
            return _RandomNumbers
        if fp is not None:
            try:
                _mod = imp.load_module('_RandomNumbers', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _RandomNumbers = swig_import_helper()
    del swig_import_helper
else:
    import _RandomNumbers
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SimpleRNG(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SimpleRNG, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SimpleRNG, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _RandomNumbers.delete_SimpleRNG
    __del__ = lambda self : None;
    def uniform(self): return _RandomNumbers.SimpleRNG_uniform(self)
SimpleRNG_swigregister = _RandomNumbers.SimpleRNG_swigregister
SimpleRNG_swigregister(SimpleRNG)

class RandomNumbersBase(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, RandomNumbersBase, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, RandomNumbersBase, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_destroy__ = _RandomNumbers.delete_RandomNumbersBase
    __del__ = lambda self : None;
RandomNumbersBase_swigregister = _RandomNumbers.RandomNumbersBase_swigregister
RandomNumbersBase_swigregister(RandomNumbersBase)

class DevRandom(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DevRandom, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DevRandom, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.DevRandom_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.DevRandom_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_DevRandom(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _RandomNumbers.delete_DevRandom
    __del__ = lambda self : None;
    def Burn(self): return _RandomNumbers.DevRandom_Burn(self)
    def UInt64(self): return _RandomNumbers.DevRandom_UInt64(self)
    def UInt32(self): return _RandomNumbers.DevRandom_UInt32(self)
    def Double(self): return _RandomNumbers.DevRandom_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.DevRandom_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.DevRandom_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.DevRandom_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.DevRandom_loadCoreState(self, *args)
    __swig_getmethods__["oneUInt64"] = lambda x: _RandomNumbers.DevRandom_oneUInt64
    if _newclass:oneUInt64 = staticmethod(_RandomNumbers.DevRandom_oneUInt64)
    __swig_getmethods__["oneUInt32"] = lambda x: _RandomNumbers.DevRandom_oneUInt32
    if _newclass:oneUInt32 = staticmethod(_RandomNumbers.DevRandom_oneUInt32)
    __swig_getmethods__["oneDouble"] = lambda x: _RandomNumbers.DevRandom_oneDouble
    if _newclass:oneDouble = staticmethod(_RandomNumbers.DevRandom_oneDouble)
DevRandom_swigregister = _RandomNumbers.DevRandom_swigregister
DevRandom_swigregister(DevRandom)

def DevRandom_defaultOptions():
  return _RandomNumbers.DevRandom_defaultOptions()
DevRandom_defaultOptions = _RandomNumbers.DevRandom_defaultOptions

def DevRandom_getCoreName():
  return _RandomNumbers.DevRandom_getCoreName()
DevRandom_getCoreName = _RandomNumbers.DevRandom_getCoreName

def DevRandom_oneUInt64():
  return _RandomNumbers.DevRandom_oneUInt64()
DevRandom_oneUInt64 = _RandomNumbers.DevRandom_oneUInt64

def DevRandom_oneUInt32():
  return _RandomNumbers.DevRandom_oneUInt32()
DevRandom_oneUInt32 = _RandomNumbers.DevRandom_oneUInt32

def DevRandom_oneDouble():
  return _RandomNumbers.DevRandom_oneDouble()
DevRandom_oneDouble = _RandomNumbers.DevRandom_oneDouble

class EmptyOptions(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, EmptyOptions, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, EmptyOptions, name)
    __repr__ = _swig_repr
    def __init__(self): 
        this = _RandomNumbers.new_EmptyOptions()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _RandomNumbers.delete_EmptyOptions
    __del__ = lambda self : None;
EmptyOptions_swigregister = _RandomNumbers.EmptyOptions_swigregister
EmptyOptions_swigregister(EmptyOptions)

class NR3Ran(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NR3Ran, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NR3Ran, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.NR3Ran_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.NR3Ran_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_NR3Ran(*args)
        try: self.this.append(this)
        except: self.this = this
    def Burn(self): return _RandomNumbers.NR3Ran_Burn(self)
    def UInt64(self): return _RandomNumbers.NR3Ran_UInt64(self)
    def UInt32(self): return _RandomNumbers.NR3Ran_UInt32(self)
    def Double(self): return _RandomNumbers.NR3Ran_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.NR3Ran_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.NR3Ran_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.NR3Ran_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.NR3Ran_loadCoreState(self, *args)
    __swig_destroy__ = _RandomNumbers.delete_NR3Ran
    __del__ = lambda self : None;
NR3Ran_swigregister = _RandomNumbers.NR3Ran_swigregister
NR3Ran_swigregister(NR3Ran)

def NR3Ran_defaultOptions():
  return _RandomNumbers.NR3Ran_defaultOptions()
NR3Ran_defaultOptions = _RandomNumbers.NR3Ran_defaultOptions

def NR3Ran_getCoreName():
  return _RandomNumbers.NR3Ran_getCoreName()
NR3Ran_getCoreName = _RandomNumbers.NR3Ran_getCoreName

class NR3Ranq1(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NR3Ranq1, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NR3Ranq1, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.NR3Ranq1_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.NR3Ranq1_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_NR3Ranq1(*args)
        try: self.this.append(this)
        except: self.this = this
    def Burn(self): return _RandomNumbers.NR3Ranq1_Burn(self)
    def UInt64(self): return _RandomNumbers.NR3Ranq1_UInt64(self)
    def UInt32(self): return _RandomNumbers.NR3Ranq1_UInt32(self)
    def Double(self): return _RandomNumbers.NR3Ranq1_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.NR3Ranq1_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.NR3Ranq1_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.NR3Ranq1_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.NR3Ranq1_loadCoreState(self, *args)
    __swig_destroy__ = _RandomNumbers.delete_NR3Ranq1
    __del__ = lambda self : None;
NR3Ranq1_swigregister = _RandomNumbers.NR3Ranq1_swigregister
NR3Ranq1_swigregister(NR3Ranq1)

def NR3Ranq1_defaultOptions():
  return _RandomNumbers.NR3Ranq1_defaultOptions()
NR3Ranq1_defaultOptions = _RandomNumbers.NR3Ranq1_defaultOptions

def NR3Ranq1_getCoreName():
  return _RandomNumbers.NR3Ranq1_getCoreName()
NR3Ranq1_getCoreName = _RandomNumbers.NR3Ranq1_getCoreName

class NR3Ranq2(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NR3Ranq2, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NR3Ranq2, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.NR3Ranq2_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.NR3Ranq2_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_NR3Ranq2(*args)
        try: self.this.append(this)
        except: self.this = this
    def Burn(self): return _RandomNumbers.NR3Ranq2_Burn(self)
    def UInt64(self): return _RandomNumbers.NR3Ranq2_UInt64(self)
    def UInt32(self): return _RandomNumbers.NR3Ranq2_UInt32(self)
    def Double(self): return _RandomNumbers.NR3Ranq2_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.NR3Ranq2_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.NR3Ranq2_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.NR3Ranq2_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.NR3Ranq2_loadCoreState(self, *args)
    __swig_destroy__ = _RandomNumbers.delete_NR3Ranq2
    __del__ = lambda self : None;
NR3Ranq2_swigregister = _RandomNumbers.NR3Ranq2_swigregister
NR3Ranq2_swigregister(NR3Ranq2)

def NR3Ranq2_defaultOptions():
  return _RandomNumbers.NR3Ranq2_defaultOptions()
NR3Ranq2_defaultOptions = _RandomNumbers.NR3Ranq2_defaultOptions

def NR3Ranq2_getCoreName():
  return _RandomNumbers.NR3Ranq2_getCoreName()
NR3Ranq2_getCoreName = _RandomNumbers.NR3Ranq2_getCoreName

class NR3Ranfib(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NR3Ranfib, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NR3Ranfib, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.NR3Ranfib_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.NR3Ranfib_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_NR3Ranfib(*args)
        try: self.this.append(this)
        except: self.this = this
    def Burn(self): return _RandomNumbers.NR3Ranfib_Burn(self)
    def Double(self): return _RandomNumbers.NR3Ranfib_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.NR3Ranfib_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.NR3Ranfib_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.NR3Ranfib_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.NR3Ranfib_loadCoreState(self, *args)
    __swig_destroy__ = _RandomNumbers.delete_NR3Ranfib
    __del__ = lambda self : None;
NR3Ranfib_swigregister = _RandomNumbers.NR3Ranfib_swigregister
NR3Ranfib_swigregister(NR3Ranfib)

def NR3Ranfib_defaultOptions():
  return _RandomNumbers.NR3Ranfib_defaultOptions()
NR3Ranfib_defaultOptions = _RandomNumbers.NR3Ranfib_defaultOptions

def NR3Ranfib_getCoreName():
  return _RandomNumbers.NR3Ranfib_getCoreName()
NR3Ranfib_getCoreName = _RandomNumbers.NR3Ranfib_getCoreName

class NR2Ran2(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, NR2Ran2, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, NR2Ran2, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.NR2Ran2_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.NR2Ran2_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_NR2Ran2(*args)
        try: self.this.append(this)
        except: self.this = this
    def Burn(self): return _RandomNumbers.NR2Ran2_Burn(self)
    def Double(self): return _RandomNumbers.NR2Ran2_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.NR2Ran2_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.NR2Ran2_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.NR2Ran2_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.NR2Ran2_loadCoreState(self, *args)
    __swig_destroy__ = _RandomNumbers.delete_NR2Ran2
    __del__ = lambda self : None;
NR2Ran2_swigregister = _RandomNumbers.NR2Ran2_swigregister
NR2Ran2_swigregister(NR2Ran2)

def NR2Ran2_defaultOptions():
  return _RandomNumbers.NR2Ran2_defaultOptions()
NR2Ran2_defaultOptions = _RandomNumbers.NR2Ran2_defaultOptions

def NR2Ran2_getCoreName():
  return _RandomNumbers.NR2Ran2_getCoreName()
NR2Ran2_getCoreName = _RandomNumbers.NR2Ran2_getCoreName

class RanLuxV32(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, RanLuxV32, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, RanLuxV32, name)
    __repr__ = _swig_repr
    __swig_getmethods__["defaultOptions"] = lambda x: _RandomNumbers.RanLuxV32_defaultOptions
    if _newclass:defaultOptions = staticmethod(_RandomNumbers.RanLuxV32_defaultOptions)
    def __init__(self, *args): 
        this = _RandomNumbers.new_RanLuxV32(*args)
        try: self.this.append(this)
        except: self.this = this
    def Burn(self): return _RandomNumbers.RanLuxV32_Burn(self)
    def Double(self): return _RandomNumbers.RanLuxV32_Double(self)
    __swig_getmethods__["getCoreName"] = lambda x: _RandomNumbers.RanLuxV32_getCoreName
    if _newclass:getCoreName = staticmethod(_RandomNumbers.RanLuxV32_getCoreName)
    def saveCoreState(self, *args): return _RandomNumbers.RanLuxV32_saveCoreState(self, *args)
    def loadCoreState(self, *args): return _RandomNumbers.RanLuxV32_loadCoreState(self, *args)
    __swig_destroy__ = _RandomNumbers.delete_RanLuxV32
    __del__ = lambda self : None;
RanLuxV32_swigregister = _RandomNumbers.RanLuxV32_swigregister
RanLuxV32_swigregister(RanLuxV32)

def RanLuxV32_defaultOptions():
  return _RandomNumbers.RanLuxV32_defaultOptions()
RanLuxV32_defaultOptions = _RandomNumbers.RanLuxV32_defaultOptions

def RanLuxV32_getCoreName():
  return _RandomNumbers.RanLuxV32_getCoreName()
RanLuxV32_getCoreName = _RandomNumbers.RanLuxV32_getCoreName

class RandomNumbers(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, RandomNumbers, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, RandomNumbers, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _RandomNumbers.new_RandomNumbers(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _RandomNumbers.delete_RandomNumbers
    __del__ = lambda self : None;
    def UInt64(self): return _RandomNumbers.RandomNumbers_UInt64(self)
    def UInt32(self): return _RandomNumbers.RandomNumbers_UInt32(self)
    def Double(self): return _RandomNumbers.RandomNumbers_Double(self)
    def Uniform(self): return _RandomNumbers.RandomNumbers_Uniform(self)
    def Exponential(self): return _RandomNumbers.RandomNumbers_Exponential(self)
    def Normal(self): return _RandomNumbers.RandomNumbers_Normal(self)
    def Gamma(self, *args): return _RandomNumbers.RandomNumbers_Gamma(self, *args)
    def GammaByMeanAndStdDev(self, *args): return _RandomNumbers.RandomNumbers_GammaByMeanAndStdDev(self, *args)
    def Poisson(self, *args): return _RandomNumbers.RandomNumbers_Poisson(self, *args)
    def Binomial(self, *args): return _RandomNumbers.RandomNumbers_Binomial(self, *args)
    def InverseCDF(self, *args): return _RandomNumbers.RandomNumbers_InverseCDF(self, *args)
    def GenerateInverseCDF(self, *args): return _RandomNumbers.RandomNumbers_GenerateInverseCDF(self, *args)
    __swig_getmethods__["defaultFilename"] = lambda x: _RandomNumbers.RandomNumbers_defaultFilename
    if _newclass:defaultFilename = staticmethod(_RandomNumbers.RandomNumbers_defaultFilename)
RandomNumbers_swigregister = _RandomNumbers.RandomNumbers_swigregister
RandomNumbers_swigregister(RandomNumbers)

def RandomNumbers_defaultFilename():
  return _RandomNumbers.RandomNumbers_defaultFilename()
RandomNumbers_defaultFilename = _RandomNumbers.RandomNumbers_defaultFilename

class SimpleRNGAdapter(SimpleRNG):
    __swig_setmethods__ = {}
    for _s in [SimpleRNG]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, SimpleRNGAdapter, name, value)
    __swig_getmethods__ = {}
    for _s in [SimpleRNG]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, SimpleRNGAdapter, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _RandomNumbers.new_SimpleRNGAdapter(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _RandomNumbers.delete_SimpleRNGAdapter
    __del__ = lambda self : None;
    def uniform(self): return _RandomNumbers.SimpleRNGAdapter_uniform(self)
SimpleRNGAdapter_swigregister = _RandomNumbers.SimpleRNGAdapter_swigregister
SimpleRNGAdapter_swigregister(SimpleRNGAdapter)

# This file is compatible with both classic and new-style classes.


