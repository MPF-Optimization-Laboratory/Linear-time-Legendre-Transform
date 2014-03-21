#include <math.h>

#include "ConjugateFunction.h"
#include "Utils.h"
#include <time.h>

float function(float *x)
{
    return x[0]*x[0] + x[1]*x[1];//abs(x[0]) + abs(x[1]);
}

int main(int argc, char* argv[])
{
    /*********************************************
     * specify the number of samples to take,
     * and how far apart they are to be spread,
     * and at what value they should start
     *********************************************/
    int   N           = 2;
    int   x_size      = 1000;
    float sample_dist = 0.05;
    float start       = -50.0f;

    /*************************************
     * Create the domain we want to run on
     *************************************/
    float *domain = new float[N * x_size];

    for (int j = 0; j < N; j++)
        constructRange(start, sample_dist, x_size, &domain[j*x_size]);

    /*************************************
     * Create the range we want to run on
     *************************************/
    int    range_size = pow((float)x_size, N);
    float *range      = new float[range_size];

    evaluateDomain(x_size, N, domain, function, range);

    /**********************************************
     * Create the slopes we would like to evaluate
     **********************************************/
    int s_size  = 2000;
    start       = -50;
    sample_dist = 0.05;

    float *slopes = new float[N * s_size];

    for (int j = 0; j < N; j++)
        constructRange(start, sample_dist, s_size, &slopes[j*s_size]);

    /**********************************************
     * Create our conjugate function and evaluate
     * at the selected points
     **********************************************/
    int    sol_size = pow((float)s_size, N);
    float *solution = new float[sol_size]; // all possible combinations of conj input

    /*********************************************
     * Initialize timers
     *********************************************/
    clock_t time_start = clock();
    double time;

    /********************************************
     * compute LLT
     ********************************************/
    linearLegendre(N, domain, range, x_size, slopes, s_size, solution);

    time = (clock() - time_start) * 1000 / CLOCKS_PER_SEC;

    /********************************************
     * print the inputs/results to a file
     ********************************************/
    //printFunction(N, x_size, domain, range, "inputs.txt");
    //printFunction(N, s_size, slopes, solution, "outputs.txt");
    printTime(time, "timing.txt", "x^2");

    /****************************************
     * Clean up allocated resources
     ****************************************/
    delete [] domain;
    delete [] range;
    delete [] slopes;
    delete [] solution;

    return 0;
}