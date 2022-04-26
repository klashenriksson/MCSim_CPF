typedef struct corr {
    double* e_buffer;
    double e_tot;
    int time_steps;
    int curr_step;
} corr_t;

corr_t corr_create(int time_steps);
void corr_step(corr_t* corr, double e);
void corr_compute(corr_t* corr, double* out_correlations, int max_correlation_steps);