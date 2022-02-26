// 2021-07-19: simkin.h is a header file for the project simkin.c

// data types
struct input_parameters
    {
    int array_size;

    double total_time;  // in seconds
    double time_step;

    int react_num;
    int comp_num;

    double tPlus;
    double tMinus;
    double T0;
    double S;
    double K1;
    double K2;
    double K3;
    double K4;
    double K5;
    };

struct variables
    {
    double *time_array;
    double **kin_arrays;
    };


// prototypes
int simulate(struct input_parameters *inpdata, struct variables *traces);
int read_input_data(FILE *inpfile, struct input_parameters *inputdata);
int allocate_memory(struct input_parameters *inputpar, struct variables *traces);
