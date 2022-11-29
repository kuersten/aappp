from distutils.core import setup, Extension
import subprocess

extra_c=subprocess.getoutput('pkg-config --cflags glib-2.0').split()
extra_l=subprocess.getoutput('pkg-config --libs glib-2.0').split()
setup(name="aappp", version="1.1", ext_modules=[Extension('aappp', ['src/aapppmodule.c', 'src/c_functions_aappp.c'], extra_link_args=extra_l, extra_compile_args=extra_c)])
