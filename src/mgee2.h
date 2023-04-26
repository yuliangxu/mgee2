/* The first half of this file is from */
/* ugee.h in the gee package */
#include <R.h>
// #include "S.h"
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R_ext/Applic.h> /* BLAS */
#include <R_ext/Linpack.h>
#include <inttypes.h>

#define PERMANENT 1
#define EPHEMERAL 0

typedef struct matrix
		{
		int nrows, ncols;
		double *data;
		int permanence;
		} MATRIX;
// class MATRIX
// {
//   int nrows, ncols;
//   double *data;
//   int permanence;
// } ;
#define ELREF( matp , s1, s2 ) ((matp)->data)+(s2)+((s1)*(matp->ncols))
#define MEL(X ,i, j) (*(ELREF( (X), (i), (j) ) ))
#define get_nelem( x ) (((x)->nrows) * ((x)->ncols))

#define malloc(n) S_alloc(n, 1)
#define calloc S_alloc
#define free(p) {p;}
#define cfree(p) {p;}
#define is_permanent( x ) (x)->permanence == PERMANENT
#define is_ephemeral( x ) (x)->permanence == EPHEMERAL
#define make_permanent( x ) (x)->permanence = PERMANENT;
#define make_ephemeral( x ) (x)->permanence = EPHEMERAL;
#define free_if_ephemeral( x ) if (is_ephemeral((x))) VC_GEE_destroy_matrix((x))

#define from_S( Sdblptr , Srowintptr , Scolintptr , Matptr ) \
Matptr = VC_GEE_create_matrix( (int)*Srowintptr, (int)*Scolintptr , EPHEMERAL ); \
{ \
int i, j, Scol, Srow; \
double *Sload; \
Scol = *Scolintptr; \
Srow = *Srowintptr; \
Sload = Sdblptr; \
for ( j = 0 ; j < Scol ; j++ ) \
	{ \
	for ( i = 0 ; i < Srow ; i++ ) \
		{ \
		MEL( Matptr , i , j ) = (double) * ( Sload ++ ); \
		} \
	} \
}

#define to_S( Matptr, Sdblptr ) \
{ \
int i, j; \
double *Sload; \
Sload = Sdblptr; \
for ( j = 0 ; j < Matptr->ncols ; j++ ) \
	{ \
	for ( i = 0 ; i < Matptr->nrows ; i++ ) \
		{ \
		* ( Sload ++ ) = MEL( Matptr , i , j ); \
		} \
	} \
}


/**************************************/
/*      Subfunctions Declarations     */
/**************************************/

static MATRIX *VC_GEE_create_matrix(),
     *VC_GEE_matcopy(),
     *VC_GEE_extract_rows(),
     *VC_GEE_matadd(),
     *VC_GEE_matsub(),
     *VC_GEE_matmult(),
     *VC_GEE_transp(),
     *VC_GEE_col_1s(),
     *VC_GEE_matabs(),
     *VC_GEE_matexp(),
     *VC_GEE_px1_times_pxq(),
     *VC_GEE_pxq_divby_px1(),
     *VC_GEE_scalar_times_matrix(),
     *VC_GEE_ident(),
     *VC_GEE_form_diag(),
     *VC_GEE_extract_cols(),
     *VC_GEE_diag_as_vec(),
     *VC_GEE_matsqrt(),
     *VC_GEE_mat1over()
     ;
static double VC_GEE_matmax(),
     VC_GEE_elsum()
     ;
static void VC_GEE_plug(),
     VC_GEE_destroy_matrix()
     ;
static int VC_GEE_split(),
     VC_GEE_nchanges()
     ;


/* The following functions are written by Z. Chen */
/* ---------------------------------------------- */

static MATRIX *get_seq1(),
     *get_rep_scalar(),
     *get_rep(),
     *get_kronecker(),
     *get_sum1row(),
     *get_sum2col(),
     *VC_GEE_matexpit(),
     *get_outer(),
     *get_rbind(),
     *get_cbind(),
     *get_cholinv(),
     *matrix_subtract(),
     *matrix_multiply(),
     *get_matrix_row()
     ;

static int get_rowindex();
static double get_max_reldif(),
              get_1_colsum(),
              get_1_rowsum();
static void get_mu_i(),
     get_lambda_i(),
     get_dmu_i_dbetaT(),
     get_dlambda_i_dbetaT(),
     get_bivar_cumuls_i(),
     get_bivar_marginals_i(),
     matrix_copyto(),
     row_replace(),
     col_replace(),
     rows_plug(),
     cols_plug(),
     matrix_addto(),
     get_matsub(),
     get_matadd(),
     get_mattransp(),
     get_matmult(),
     cholinv(),
     add_outer_colvec_to(),
     scalar_times_matrix(),
     matrix_elem_mult(),
     matrix_row_mult(),
     matrix_col_mult(),
     get_estfun(),
     get_dvd(),
     fisherscoring(),
     set_zero(),
     get_dpXt_i_dvpT_l()
     ;

static MATRIX *get_dYtil_i_dgT(),
     *get_dZtil_i_dgT(),
     *form_matrix();

void Cmgee2(),
     Cgetmgee2_i(),
     Cgetordgee2_i()
     ;

