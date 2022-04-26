#include "corr.h"

#include <stdio.h>
#include <stdlib.h>

corr_t corr_create(int time_steps)
{
    corr_t r;
    r.e_buffer = malloc(sizeof(double)*time_steps);
    r.time_steps = time_steps;
    r.curr_step = 0;
    r.e_tot = 0;
    return r;
}

void corr_step(corr_t* corr, double e)
{
    corr->e_buffer[corr->curr_step++] = e;
    corr->e_tot += e;
}

void corr_compute(corr_t* corr, double* out_correlations, int max_correlation_steps)
{
    double e_avg = corr->e_tot / corr->time_steps;
    for (int i = 0; i < max_correlation_steps; i++)
    {
        double C = 0.;
        for (int data_i = 0; data_i < corr->time_steps-i; data_i++)
        {
            double e_0 = corr->e_buffer[data_i];
            double e_t = corr->e_buffer[data_i+i];
            C += (e_0-e_avg)*(e_t-e_avg);
        }

        C /= (corr->time_steps-i);
        out_correlations[i] = C;
    }
}