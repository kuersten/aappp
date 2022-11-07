#!bash

#CFLAGS="-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include" python3 setup.py sdist
python3 setup.py sdist	#compile package
pip3 install ./dist/aappp-1.0.tar.gz	#install package
