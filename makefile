CC = mpic++
CFLAGS = -Wall -std=c++11  

SRCS = main.cpp grad_conj.cpp produit.cpp fonction.cpp charge.cpp

OBJS = $(SRCS:.cpp=.o)

TARGET = run

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

launch: $(TARGET)
	mpirun -np 4 ./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) *.dat

help:
	@echo "Usage:"
	@echo "  make       : Compile l'exÃ©cutable $(TARGET)"
	@echo "  make clean : Supprime les fichiers objets et l'exÃ©cutable"
