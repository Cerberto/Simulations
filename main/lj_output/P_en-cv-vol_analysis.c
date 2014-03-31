#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main (int argc, char *argv[]) {
	if (argc<2) {
		printf("Error : missing output_file and number of lines.\n\n");
		exit(EXIT_FAILURE);
	}

	double press, en, en1, cv, cv1, vol, vol1, spare;
	double en_err, cv_err, vol_err;
	int counter, lines;

	char *output_filename;
		output_filename = malloc(100*sizeof(char));
	FILE *output;
	output = fopen(output_filename, "a");

	/*
		Read the input file and compute the averages to write on output
	*/
	en  = 0;
	cv  = 0;
	vol = 0;
	counter = 0;
	while(counter < lines) {
		scanf("%lf",&press);
		scanf("%lf",&en1);
		scanf("%lf",&spare);
		scanf("%lf",&cv1);
		scanf("%lf",&spare);
		scanf("%lf",&vol1);
		scanf("%lf",&spare);

		en  += en1;
		cv  += cv1;
		vol += vol1;
		en_err  += en1*en1;
		cv_err  += cv1*cv1;
		vol_err += vol1*vol1;
		counter++;
	}

	en  /= (double)counter;
	cv  /= (double)counter;
	vol /= (double)counter;
	en_err  /= (double)counter;
	cv_err  /= (double)counter;
	vol_err /= (double)counter;
	en_err  += -en*en;
	cv_err	+= -cv*cv;
	vol_err += -vol*vol;
	en_err  /= (double)(counter - 1);
	cv_err  /= (double)(counter - 1);
	vol_err /= (double)(counter - 1);

	en_err  = sqrt(en_err);
	cv_err  = sqrt(cv_err);
	vol_err = sqrt(vol_err);

	/* ********************* */

	fprintf(output, "%.3lf\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\t%.7e\n", \
		press, en, en_err, cv, cv_err, vol, vol_err);

	exit(EXIT_SUCCESS);
}
