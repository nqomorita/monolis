#ifndef MONOLIS_WRAPPER_ML_H
#define MONOLIS_WRAPPER_ML_H

extern void monolis_ml_get_nlocal_(int *id, int *nlocal, int *nlocal_allcolumns, int *ierr);
extern void monolis_ml_get_coord_(int *id, double x[], double y[], double z[], int *ierr);
extern void monolis_ml_get_rbm_(int *id, double rbm[], int *ierr);
extern void monolis_ml_get_loglevel_(int *id, int *level);
extern void monolis_ml_get_opt_(int *id, int opt[], int *ierr);
extern void monolis_ml_set_opt_(int *id, int opt[], int *ierr);

#endif
