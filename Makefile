SOURCES := 	./src/utils.c \
		./src/domain.c \
		./src/phys_const.c \
		./src/cosmology.c \
		./src/correlation_function.c \
		./src/cross_correlation_function.c \
		./src/main.c

OBJECTS := $(SOURCES:.c=.o)
DOBJECTS := $(SOURCES:.c=.d)
EXECUTABLE := pcorrfunc

USE-MPI=YES
# USE-SPRNG=YES
# USE-DEBUG-CORRFUNC=YES
# BUILD-LIB=YES

include common.mk


.PHONY: all clean clena celan celna

ifdef BUILD-LIB
all: $(SOURCES) $(EXECUTABLE)
	$(CC) -shared -o libpcf.so $(OBJECTS)
else
all: $(SOURCES) $(EXECUTABLE)
endif

celan celna clena:clean


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf $(OBJECTS) $(DOBJECTS) $(EXECUTABLE)
