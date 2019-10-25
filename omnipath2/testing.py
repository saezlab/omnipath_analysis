import pypath
from pypath.main import PyPath

# Setting up working environment
cachedir = '/home/nico/pypath_cache'

if os.getcwd().endswith('omnipath2'):
    os.chdir('omnipath2')

if not os.path.exists(cachedir):
    os.makedirs(cachedir)

pypath.settings.setup(cachedir=cachedir)

pa = PyPath()

#pa.init_network()
pa.init_network(pfile=os.path.join(cachedir, 'network.pickle'))

#pa.save_network(pfile=os.path.join(cachedir, 'network.pickle'))

# PPI load methods

load_methods = [m for m in dir(pa) if m.startswith('load')]
load_methods.remove('load_ielm')

len(load_methods)

ex_count = 0

for meth in load_methods:
    print('\n' + '=' * 70)
    print('Calling load method %s' % meth)

    try:
        getattr(pa, meth)()
        print('%s loaded successfully (?)' % meth)

    except Exception as err:
        ex_count += 1
        print(err)

ex_count/len(load_methods)

# dataio module

totest = [m for m in dir(dataio)
          if (str(type(getattr(dataio, m))) == "<class 'function'>"
              and not m.startswith('_'))]

len(totest)

ex_count = 0

for meth in totest:
    print('\n' + '=' * 70)
    print('Calling load method %s' % meth)

    try:
        getattr(dataio, meth)()
        print('%s loaded successfully (?)' % meth)

    except Exception as err:
        ex_count += 1
        print(err)

ex_count/len(totest)
