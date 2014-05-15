

parcel_aux.so: parcel_aux.pyx
	python setup.py build_ext --inplace
	
clean:
	rm *.pyc *.so *.c