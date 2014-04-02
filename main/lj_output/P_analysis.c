#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[]) {
	if (argc<2) {
		printf("Error : missing number of lines of input files and input files.\n\n");
		exit(EXIT_FAILURE);
	}

	double press, en, en1, cv, cv1, vol, vol1, spare;
	double en_err, cv_err, vol_err;
	int counter, lines, ifile;

	char *input_filename;
		input_filename  = malloc(100*sizeof(char));
	FILE *input, *output;
		output = fopen("P_many-sims", "a");

	lines = atoi(argv[1]);

	for(ifile=2; ifile<argc; ifile++) {

		input_filename = argv[ifile];
		input = fopen(input_filename, "r");

		/*
	 	*	Read the input file and compute the averages to write on output
		**/
		en  = 0;
		cv  = 0;
		vol = 0;
		counter = 0;
		while(counter < lines) {
			fscanf(input,"%lf",&press);
			fscanf(input,"%lf",&en1);
			fscanf(input,"%lf",&spare);
			fscanf(input,"%lf",&cv1);
			fscanf(input,"%lf",&spare);
			fscanf(input,"%lf",&vol1);
			fscanf(input,"%lf",&spare);

			en  += en1;
			cv  += cv1;
			vol += vol1;
			en_err  += en1*en1;
			cv_err  += cv1*cv1;
			vol_err += vol1*vol1;
			counter++;
		}

		en  /= (double)lines;
		cv  /= (double)lines;
		vol /= (double)lines;
		en_err  /= (double)lines;
		cv_err  /= (double)lines;
		vol_err /= (double)lines;
		en_err  += -en*en;
		cv_err	+= -cv*cv;
		vol_err += -vol*vol;
		en_err  /= (double)(lines - 1);
		cv_err  /= (double)(lines - 1);
		vol_err /= (double)(lines - 1);

		en_err  = sqrt(en_err);
		cv_err  = sqrt(cv_err);
		vol_err = sqrt(vol_err);

		/* ********************* */

		fprintf(output, "%.3lf\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", \
			press, en, en_err, cv, cv_err, vol, vol_err);

		fclose(input);
	}

	exit(EXIT_SUCCESS);
}
