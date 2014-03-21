#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <string.h>

/***************************************************
 * allows you to construct a set of x values
 *
 * input:
 * start: the initial value to be included in the range
 * incr:  the size of each step
 * count: the number of elements to be added
 *
 * in/out:
 * data: a pointer to the location in memory to
 *       construct the set
 ***************************************************/
void constructRange(float start, float incr, int count, float *data)
{
    for (int i = 0; i < count; i++)
    {
        data[i] = start;
        start += incr;
    }
}

/**************************************************************
 * allows you to construct the matrix of function values on
 * the combinations of the domain
 *
 * THIS IS NOT EFFICIENT, it would be much better to pass
 * the function into the conjugate function once implemented
 *
 * input:
 * size:   the number of elements per x_i
 * N:      the multiplicity of the function (R^N)
 * domain: a pointer to the array of x values formatted
 *         x1[0], x1[1],... x2[0], x2[1],... xn[0], xn[1]
 * func:   a pointer to the function that will evaluate
 *         a set of x inputs
 *
 * output:
 * store: the location in memory to store the output
 ****************************************************************/
void evaluateDomain(int size, int N, float *domain, float (*func)(float *), float *store)
{
    int     y_size = pow((float)size, N);
    float **value  = new float*[N]; // pointers to the start of each of the respective x_i vectors
    int    *index  = new int[N];    // the current offsets of the x values being evaluated
    float  *X      = new float[N];  // the conglomerated values at which to evaluate the function

    // get the pointers to each of the parameter arrays
    for (int i = 0; i < N; i++)
        value[i] = &domain[i*size];

    // set the index array to zero
    memset(index, 0, sizeof(int)*N);

    for (int i = 0; i < y_size; i++)
    {
        for (int k = 0; k < N; k++)
            X[k] = value[k][index[k]];
        store[i] = func(X);

        // get the next combination of x values f(x1, ..., x_n)
        for (int j = N-1; j >= 0 && ++index[j] == size; j--)
            index[j] = 0;
    }

    delete [] index;
    delete [] X;
}

void printTime(double time, char *filename, char *funcname = NULL, bool append=true)
{
    FILE *f = NULL;
    
    if (append)
        f = fopen(filename, "a");
    else
        f = fopen(filename,"w+");

    fprintf(f, "Func %s, time %lf\n", funcname, time);
}

/*************************************************************
 * input:
 * N:      the multiplicity of the function (R^N)
 * size:   the number of components per vector
 * domain: a pointer to the array of x values formatted
 *         x1[0], x1[1],... x2[0], x2[1],... xn[0], xn[1]
 * func:   a pointer to the function that will evaluate
 *         a set of x inputs
 * filename: the name of the file to write out to
 *************************************************************/
void printFunction(int N, int x_size, float *domain, float *range, char *filename)
{
    FILE *f = fopen(filename,"w+");

    int     y_size = pow((float)x_size, N);
    float **value  = new float*[N]; // pointers to the start of each of the respective x_i vectors
    int    *index  = new int[N];    // the current offsets of the x values being evaluated
    float  *X      = new float[N];  // the conglomerated values at which to evaluate the function

    // get the pointers to each of the parameter arrays
    for (int i = 0; i < N; i++)
        value[i] = &domain[i*x_size];

    // set the index array to zero
    memset(index, 0, sizeof(int)*N);

    for (int i = 0; i < y_size; i++)
    {
        int incr = 0;
        for (int k = 0; k < N; k++)
        {
            X[k] = value[k][index[k]];
            fprintf(f, "%f ", X[k]);
            incr += index[k] * pow((double)x_size, (N-1 -k));
        }
        fprintf(f, "%f\n", range[incr]);
        fflush(f);
        // get the next combination of x values f(x1, ..., x_n)
        for (int j = 0; j < N && ++index[j] == x_size; j++)
        {
            fprintf(f, "\n");
            index[j] = 0;
        }
    }

    delete [] index;
    delete [] X;
    fclose(f);
}

#endif //UTILS_H