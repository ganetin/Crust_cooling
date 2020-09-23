NAME = runeos
F77 = gfortran -O 
#F77 = frt -O 
EXEC = $(NAME)

OBJS =  main.o eos.o newt.o eos_tye_np.o 

$(NAME): $(OBJS)
	@echo building executable
	$(F77) -o $(EXEC) $(OBJS)

.f.o:
	@echo building $*.o
	$(F77) -c $<

clean:	
	rm -f *.o runeos *~
