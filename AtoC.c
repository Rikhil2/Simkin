#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mem.h>
//test this
#include <assert.h>

#include "AtoC.h"
#include "algebra_d.h"


struct input_parameters inputdata;
struct variables ttraces;

int main(void)
{


    FILE *input_file;
    //can't input data
    //assert((!read_input_data(input_file, &inputdata)));
    //two possible problems;

    if (!read_input_data(input_file, &inputdata))
        {
        puts("input file not found\n");
        puts("press any key\n");
        getchar();
        exit(0);
        };

    if (!allocate_memory(&inputdata, &ttraces))
        {
        puts("memory allocation error\n");
        puts("press any key\n");
        getchar();
        exit(0);
        };

    printf("Memory allocated successfully: %d arrays %d bytes each\n\n", inputdata.comp_num, inputdata.array_size*sizeof(double));

    simulate(&inputdata, &ttraces);


    vectprint(ttraces.time_array, ttraces.kin_arrays[0], inputdata.array_size, "S.dat");
    vectprint(ttraces.time_array, ttraces.kin_arrays[1], inputdata.array_size, "T plus.dat");
    vectprint(ttraces.time_array, ttraces.kin_arrays[2], inputdata.array_size, "T minus.dat");
    vectprint(ttraces.time_array, ttraces.kin_arrays[3], inputdata.array_size, "T zero.dat");

    getchar();
    return 0;
}


// reading input parameters
int read_input_data(FILE *inpfile, struct input_parameters *inpdata)
{
char buffer[1000];

if (!(inpfile=fopen("AtoC.par", "rt"))) return(0);

fgets(buffer, 500, inpfile);// empty line
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->total_time)); //scans the first few characters
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->time_step));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%d", &(inpdata->react_num));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%d", &(inpdata->comp_num));
fgets(buffer, 500, inpfile);


fgets(buffer, 500, inpfile);
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->S));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->T0));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->tPlus));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->tMinus));
fgets(buffer, 500, inpfile);

fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->K1));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->K2));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->K3));
fgets(buffer, 500, inpfile);

sscanf(buffer, "%lf", &(inpdata->K4));


/*
fgets(buffer, 500, inpfile);
sscanf(buffer, "%lf", &(inpdata->K5));
*/
/*
inpdata->K1=vc(inpdata->comp_num);
inpdata->K2=vc(inpdata->comp_num);
inpdata->K3=vc(inpdata->comp_num);
inpdata->K4=vc(inpdata->comp_num);
inpdata->K5=vc(inpdata->comp_num);
*/
//just in case; I did it manually
/*
for (i=0; i<inpdata->comp_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->S[i]));
    }

for (i=0; i<inpdata->comp_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->tPlus[i]));
    }
for (i=0; i<inpdata->comp_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->tMinus[i]));
    }
for (i=0; i<inpdata->comp_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->T0[i]));
    }

fgets(buffer, 500, inpfile);

for (i=0; i<inpdata->react_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->K1[i]));
    }
for (i=0; i<inpdata->react_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->K2[i]));
    }
for (i=0; i<inpdata->react_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->K3[i]));
    }
for (i=0; i<inpdata->react_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->K4[i]));
    }
for (i=0; i<inpdata->react_num; ++i)
    {
    fgets(buffer, 500, inpfile);
    sscanf(buffer, "%lf", &(inpdata->K5[i]));
    }
*/
fclose(inpfile);

return(1);
}


// memmory allocation
int allocate_memory(struct input_parameters *inputpar, struct variables *traces)
{

inputpar->array_size = 1 + inputpar->total_time/inputpar->time_step;

traces->time_array=vc(inputpar->array_size);
traces->kin_arrays=mc(inputpar->array_size, inputpar->comp_num);

return(1);
}


// performs actual calculations
int simulate(struct input_parameters *inpdata, struct variables *traces)
{
int i;
// starting conditions
traces->time_array[0]=0.0;

traces->kin_arrays[0][0]=inpdata->S; //each gets its own array, and starts at 0, so it can be changed
traces->kin_arrays[1][0]=inpdata->tPlus;
traces->kin_arrays[2][0]=inpdata->tMinus;
traces->kin_arrays[3][0]=inpdata->T0;
//problem with initialization of Ks
traces->kin_arrays[4][0]=inpdata->K1;//rate constants are constant; they do not change
traces->kin_arrays[5][0]=inpdata->K2;
traces->kin_arrays[6][0]=inpdata->K3;
traces->kin_arrays[7][0]=inpdata->K4;
traces->kin_arrays[8][0]=inpdata->K5;

//not sure what this does
//traces->kin_arrays[4][0]=traces->kin_arrays[2][0]=0.0;

for (i=1; i<inpdata->array_size; ++i)
    {
    traces->time_array[i]=inpdata->time_step*i;

    traces->kin_arrays[0][i] = traces->kin_arrays[0][i-1] + (((traces->kin_arrays[4][0]*traces->kin_arrays[3][i]) + (traces->kin_arrays[5][0]*traces->kin_arrays[1][i]) + (traces->kin_arrays[5][0]*traces->kin_arrays[2][i]))*inpdata->time_step);
    traces->kin_arrays[1][i] = traces->kin_arrays[1][i-1] + (((traces->kin_arrays[5][0]*traces->kin_arrays[0][i]) + (traces->kin_arrays[5][0]*traces->kin_arrays[3][i]) -2*(traces->kin_arrays[5][0]*traces->kin_arrays[1][i]))*inpdata->time_step);
    traces->kin_arrays[2][i] = traces->kin_arrays[2][i-1] + (((traces->kin_arrays[5][0]*traces->kin_arrays[0][i]) + (traces->kin_arrays[5][0]*traces->kin_arrays[3][i]) -2*(traces->kin_arrays[5][0]*traces->kin_arrays[2][i]))*inpdata->time_step);
    traces->kin_arrays[3][i] = traces->kin_arrays[3][i-1] + (((traces->kin_arrays[4][0]*traces->kin_arrays[0][i]) - (2*traces->kin_arrays[5][0] + traces->kin_arrays[4][0])*traces->kin_arrays[3][i] + (traces->kin_arrays[5][0]*traces->kin_arrays[1][i]) + (traces->kin_arrays[5][0]*traces->kin_arrays[2][i]))*inpdata->time_step);
    }

return(1);
}
