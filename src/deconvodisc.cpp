// Poisson (multinomial) data from known convolution of smooth discrete distro
#define TMB_LIB_INIT R_init_deconvodisc
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // IMO this should be IMATRIX, but...
  DATA_MATRIX( count);         // Observations (one dataset per column)
  DATA_MATRIX( Desmat_fix);  // Random effect design matrix
  DATA_MATRIX( Desmat_rand);  // Fixed effect design matrix
  DATA_MATRIX( chol_S);   // full-rank bit of S(moother)
  DATA_MATRIX( M);        // convolution
  PARAMETER_MATRIX(u);    // Random effects vectors (one dataset per col)
  PARAMETER_MATRIX(beta); // Fixed effects vectors
  PARAMETER(logsdu);      // Random effect standard deviations

  // Distribution of random effect (u):
  Type lglk = 0;

  matrix<Type> u_indept = chol_S * u;
  lglk += dnorm( u_indept.vec(), Type(0), exp(logsdu), true).sum();

  // Distribution of obs given random effects (x|u):
  // matrix<Type> E_unconv = exp( Desmat_fix * beta + Desmat_rand * u); doesn't work :/
  // Two steps needed; but that's not all...
  matrix<Type> log_E_unconv = Desmat_fix * beta + Desmat_rand * u;
  matrix<Type> E_unconv = exp( log_E_unconv.array()); // sans array() it f***s up dims :/

/*
  vector<int> d_E_unconv(2);
  d_E_unconv(0) = E_unconv. rows();
  d_E_unconv(1) = E_unconv. cols();
  REPORT( d_E_unconv);
*/

  matrix<Type> E_count = M * E_unconv;
  
  // Probabilities (the real goal here)
  matrix<Type> E_norm( E_unconv.rows(), E_unconv.cols());
  for( int j=0; j < E_unconv.cols(); j++) {
    Type inv_jsum = Type( 1.0) / E_unconv. col( j). sum();
    
    for( int i=0; i < E_unconv.rows(); i++) {
      E_norm( i, j) = E_unconv( i, j) * inv_jsum;
    };
  };
  ADREPORT( E_norm);

  lglk += dpois( count.vec(), E_count.vec(), true).sum();

  return -lglk;
}

