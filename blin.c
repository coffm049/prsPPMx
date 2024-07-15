#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <time.h>


void gibbs_sampler(const gsl_vector *y, const gsl_matrix *X, const gsl_vector *beta0, double alpha0, double beta0_prior, int n_iter, int burn_in, gsl_matrix *beta_samples, gsl_vector *sigma2_samples) {
    int n = y->size;
    int p = X->size2;
    
    // Initialize V0 to 100 * identity matrix
    gsl_matrix *V0 = gsl_matrix_alloc(p, p);
    gsl_matrix_set_identity(V0);
    gsl_matrix_scale(V0, 100.0);
    
    // Initialize storage matrices
    gsl_vector *beta = gsl_vector_alloc(p);
    gsl_vector_memcpy(beta, beta0);
    double sigma2 = 1.0;
    
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(NULL));
    
    // Gibbs sampling
    for (int i = 0; i < n_iter; i++) {
        // Sample beta | y, sigma2
        gsl_matrix *XtX = gsl_matrix_alloc(p, p);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, X, X, 0.0, XtX);
        gsl_matrix *V0_inv = gsl_matrix_alloc(p, p);
        gsl_matrix_memcpy(V0_inv, V0);
        gsl_linalg_cholesky_decomp(V0_inv);
        gsl_linalg_cholesky_invert(V0_inv);
        gsl_matrix_add(XtX, V0_inv);
        
        gsl_matrix *Vn = gsl_matrix_alloc(p, p);
        gsl_matrix_memcpy(Vn, XtX);
        gsl_linalg_cholesky_decomp(Vn);
        gsl_linalg_cholesky_invert(Vn);
        
        gsl_vector *Xty = gsl_vector_alloc(p);
        gsl_blas_dgemv(CblasTrans, 1.0, X, y, 0.0, Xty);
        gsl_vector *V0_inv_beta0 = gsl_vector_alloc(p);
        gsl_blas_dgemv(CblasNoTrans, 1.0, V0_inv, beta0, 0.0, V0_inv_beta0);
        gsl_vector_add(Xty, V0_inv_beta0);
        
        gsl_vector *betan = gsl_vector_alloc(p);
        gsl_blas_dgemv(CblasNoTrans, 1.0, Vn, Xty, 0.0, betan);
        
        gsl_matrix *L = gsl_matrix_alloc(p, p);
        gsl_matrix_memcpy(L, Vn);
        gsl_linalg_cholesky_decomp(L);
        
        for (int j = 0; j < p; j++) {
            gsl_vector_set(beta, j, gsl_ran_gaussian(r, sqrt(sigma2 * gsl_matrix_get(Vn, j, j))) + gsl_vector_get(betan, j));
        }
        
        // Sample sigma2 | y, beta
        double alpha_n = alpha0 + n / 2.0;
        gsl_vector *residuals = gsl_vector_alloc(n);
        gsl_vector_memcpy(residuals, y);
        gsl_blas_dgemv(CblasNoTrans, -1.0, X, beta, 1.0, residuals);
        
        double ssq = 0.0;
        for (int k = 0; k < n; k++) {
            ssq += gsl_vector_get(residuals, k) * gsl_vector_get(residuals, k);
        }
        
        double beta_n = beta0_prior + 0.5 * ssq;
        sigma2 = 1.0 / gsl_ran_gamma(r, alpha_n, 1.0 / beta_n);
        
        // Store samples
        if (i >= burn_in) {
            for (int j = 0; j < p; j++) {
                gsl_matrix_set(beta_samples, i - burn_in, j, gsl_vector_get(beta, j));
            }
            gsl_vector_set(sigma2_samples, i - burn_in, sigma2);
        }
        
        // Clean up
        gsl_matrix_free(XtX);
        gsl_matrix_free(V0_inv);
        gsl_matrix_free(Vn);
        gsl_vector_free(Xty);
        gsl_vector_free(V0_inv_beta0);
        gsl_vector_free(betan);
        gsl_matrix_free(L);
        gsl_vector_free(residuals);
    }
    
    gsl_matrix_free(V0);
    gsl_vector_free(beta);
    gsl_rng_free(r);
}

int main() {
    // Example usage
    int n = 500; // number of observations
    int p = 2; // number of predictors
    
    gsl_vector *y = gsl_vector_alloc(n);
    gsl_matrix *X = gsl_matrix_alloc(n, p);
    
    // Simulate data
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(NULL));
    
    double true_beta[] = {2.0, -1.0};
    double true_sigma2 = 4.0;
    
    for (int i = 0; i < n; i++) {
        double x1 = gsl_ran_gaussian(r, 1.0);
        double x2 = gsl_ran_gaussian(r, 1.0);
        double epsilon = gsl_ran_gaussian(r, sqrt(true_sigma2));
        
        gsl_matrix_set(X, i, 0, x1);
        gsl_matrix_set(X, i, 1, x2);
        gsl_vector_set(y, i, true_beta[0] * x1 + true_beta[1] * x2 + epsilon);
    }
    
    gsl_vector *beta0 = gsl_vector_alloc(p);
    for (int j = 0; j < p; j++) {
        gsl_vector_set(beta0, j, 0.0);
    }
    
    double alpha0 = 2.0;
    double beta0_prior = 1.0;
    int n_iter = 1000;
    int burn_in = 500;
    
    gsl_matrix *beta_samples = gsl_matrix_alloc(n_iter - burn_in, p);
    gsl_vector *sigma2_samples = gsl_vector_alloc(n_iter - burn_in);
    
    gibbs_sampler(y, X, beta0, alpha0, beta0_prior, n_iter, burn_in, beta_samples, sigma2_samples);
    
    // Output samples
    for (int i = 0; i < n_iter - burn_in; i++) {
        printf("Iteration %d: beta = [", i + 1);
        for (int j = 0; j < p; j++) {
            printf("%f", gsl_matrix_get(beta_samples, i, j));
            if (j < p - 1) {
                printf(", ");
            }
        }
        printf("], sigma2 = %f\n", gsl_vector_get(sigma2_samples, i));
    }
    
    // Clean up
    gsl_vector_free(y);
    gsl_matrix_free(X);
    gsl_vector_free(beta0);
    gsl_matrix_free(beta_samples);
    gsl_vector_free(sigma2_samples);
    gsl_rng_free(r);
    
    return 0;
}

