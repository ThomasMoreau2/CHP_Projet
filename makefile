# Définir le compilateur et les options
CC = g++
CFLAGS = -Wall -std=c++11  # Utilisez cette option pour le standard C++11

# Fichiers sources
SRCS = main.cpp grad_conj.cpp produit.cpp fonction.cpp

# Fichiers objets (remplacer .cpp par .o)
OBJS = $(SRCS:.cpp=.o)

# Nom de l'exécutable
TARGET = run

# Règle par défaut (tout)
all: $(TARGET)

# Règle pour lier l'exécutable
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $(TARGET)

# Règle pour compiler chaque fichier source en fichier objet
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Règle pour nettoyer les fichiers générés (objets et exécutable)
clean:
	rm -f $(OBJS) $(TARGET)

# Règle pour afficher les informations d'aide
help:
	@echo "Usage:"
	@echo "  make       : Compile l'exécutable $(TARGET)"
	@echo "  make clean : Supprime les fichiers objets et l'exécutable"
