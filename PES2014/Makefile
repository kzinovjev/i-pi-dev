################################################################################
#
#   This Makefile will compile the potential energy surface module to expose it
#   to Python using f2py.
#
################################################################################

F2PY=f2py

.PHONY: all build pyf clean

all: build

build:
	$(F2PY) -c ch4oh2014.f utility.f ch4oh2014.pyf

pyf:
	$(F2PY) -m PES -h ch4oh2014.pyf ch4oh2014.f --overwrite-signature

clean:
	rm -f *.o *.so
