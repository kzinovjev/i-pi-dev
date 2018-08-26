include $(AMBERHOME)/config.h

.PHONY: clean

sander_driver: ../sockets.o ../fsockets.o sander_driver.o
	$(FC) $(FFLAGS) -I$(INCDIR) -L$(LIBDIR) $(LDFLAGS) -o $@ $^ -lsander

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod sander_driver