#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[]) {
	
	/* Errore con meno di 3 argomenti (oltre al nome del programma) */
	if (argc < 4) {
		printf("Error : missing arguments!\n");
		exit(EXIT_FAILURE);
	}

	double press, en, en_err, cv, cv_bis, cv_err, vol, vol_err, spare;

	/* Primo argomento è un file di input dove leggo energia, cv e volume */
	char *input1_file;
	input1_file = malloc(100*sizeof(char));
	input1_file = argv[1];
	FILE *i1 = fopen(input1_file, "r");

	/* Secondo argomento è un file di input dove leggo solo cv (da sostituire a quello del primo input) */
	char *input2_file;
	input2_file = malloc(100*sizeof(char));
	input2_file = argv[2];
	FILE *i2 = fopen(input2_file, "r");

	/* Terzo argomento è il numero di righe da leggere */
	int lines;
	lines = atoi(argv[3]);

	/* Apro un file di output */
	char *output;
	output = malloc(100*sizeof(char));
	FILE *out;

	/* Leggo  e stampo */
	int counter = 0;
	fscanf(i1,"%lf",&press);
	sprintf(output, "P_en-cv-vol_ok_%.3lf", press);
	out = fopen(output, "w");
	while(counter < lines) {
		if(counter>0)
			fscanf(i1,"%lf",&spare);
		fscanf(i1,"%lf",&en);
		fscanf(i1,"%lf",&en_err);
		fscanf(i1,"%lf",&spare);
		fscanf(i1,"%lf",&spare);
		fscanf(i1,"%lf",&vol);
		fscanf(i1,"%lf",&vol_err);

		fscanf(i2,"%lf",&press);
		fscanf(i2,"%lf",&cv_bis); 
		fscanf(i2,"%lf",&cv_err);

		fprintf(out, "%.3lf\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", \
			press, en, en_err, cv_bis, cv_err, vol, vol_err);

		counter++;
	}

	exit(EXIT_SUCCESS);
}