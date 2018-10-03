SOURCES := 	./utils.c \
		./domain.c \
		./phys_const.c \
		./cosmology.c \
		./correlation_function.c \
		./cross_correlation_function.c \
		./main.c

OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := pcorrfunc

USE-MPI=YES
# USE-SPRNG=YES
# USE-DEBUG-CORRFUNC=YES
 
include common.mk


.PHONY: all clean clena celan celna

all: $(SOURCES) $(EXECUTABLE)

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
