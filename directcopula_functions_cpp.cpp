#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <string>
#include <sstream>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace arma;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double sum_na_vec(vec X) {
    double sum = 0;
    for (int i = 0; i < X.n_elem; i++) {
            sum += X(i);
    }
    return sum;
}

// [[Rcpp::export]]
double prod_na_vec(vec X) {
  double sum = 1;
  for (int i = 0; i < X.n_elem; i++) {
    sum = sum * X(i);
  }
  return sum;
}

// [[Rcpp::export]]
vec row_sums(mat beta){
  int mat_nrow=beta.n_rows;
  vec ret = zeros<vec>(mat_nrow);
  for(int i=0;i<mat_nrow;i++){
    ret(i)=sum_na_vec(beta.row(i).t());
  }
  return ret;
}


// [[Rcpp::export]]
vec col_sums(mat beta){
  int mat_ncol=beta.n_cols;
  vec ret = zeros<vec>(mat_ncol);
  for(int i=0;i<mat_ncol;i++){
    ret(i)=sum_na_vec(beta.col(i));
  }
  return ret;
}


// [[Rcpp::export]]
vec col_prods(mat beta){
  int mat_ncol=beta.n_cols;
  vec ret = zeros<vec>(mat_ncol);
  for(int i=0;i<mat_ncol;i++){
    ret(i)=prod_na_vec(beta.col(i));
  }
  return ret;
}

// [[Rcpp::export]]
vec row_prods(mat beta){
  int mat_nrow=beta.n_rows;
  vec ret = zeros<vec>(mat_nrow);
  for(int i=0;i<mat_nrow;i++){
    ret(i)=prod_na_vec(beta.row(i).t());
  }
  return ret;
}

// [[Rcpp::export]]
mat row_prods_mat(mat beta){
  int mat_nrow=beta.n_rows;
  mat ret = zeros<mat>(mat_nrow,1);
  for(int i=0;i<mat_nrow;i++){
    ret(i,0)=prod_na_vec(beta.row(i).t());
  }
  return ret;
}


// [[Rcpp::export]]
mat mylogit_cpp(mat x){
      mat y=1/(1+exp(-x));
      return y;
}

// [[Rcpp::export]]
mat mylogitpdf_cpp(mat x){
      mat y=exp(-x)/(pow( 1 + exp(-x),2));
      return y;
}

// [[Rcpp::export]]
mat gen_sine_cpp(mat u){
      mat ret = sin(datum::pi * u)/datum::pi;
      return ret;
}

// [[Rcpp::export]]
mat dgen_sine_cpp(mat u){
      mat ret = cos(datum::pi * u);
      return ret;
}

// [[Rcpp::export]]
mat gen_power_cpp(mat u, double n){
      mat ret = 1-pow(pow(u,n)+pow(1-u,n),1/n);
      return ret;
}

// [[Rcpp::export]]
mat dgen_power_cpp(mat u, double n){
      mat ret = (pow(1-u,n-1)-pow(u,n-1)) % pow(pow(u,n)+pow(1-u,n),1/n-1);
      return ret;
}


// [[Rcpp::export]]
mat el_pow(mat x,mat p){
    int mat_ncol=x.n_cols;
    int mat_nrow=x.n_rows;
    mat ret = zeros<mat>(mat_nrow,mat_ncol);
    for(int i=0;i<mat_nrow;i++){
            for(int j=0;j<mat_ncol;j++){
                    ret(i,j) = pow(x(i,j),p(i,j));
            }
    }
    return ret;
}

// [[Rcpp::export]]
mat likelihood_estM_power_cpp(mat x,mat pur,mat attri,mat price,vec beta,vec theta,mat tau_mat,int nalt,int nbeta,int ntheta,double n){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 /(x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat ggg = el_pow(gen_power_cpp(mylogit_cpp(G),n)/mylogit_cpp(G),1-pur) % el_pow(dgen_power_cpp(mylogit_cpp(G),n),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(mylogit_cpp(G),1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(mylogitpdf_cpp(G),pur));
  return ret;

}


// [[Rcpp::export]]
mat likelihood_estM_power_cpp_norm(mat x,mat pur,mat attri,mat price,vec beta,vec theta,mat tau_mat,int nalt,int nbeta,int ntheta,double n){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 /(x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat G_dn = mat(nobs,nalt);
  mat G_pn = mat(nobs,nalt);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G_dn(i,j) = R::dnorm(G(i,j),0.0,1.0,false);
      G_pn(i,j) = R::pnorm(G(i,j),0.0,1.0,true,false);
    }
  }

  mat ggg = el_pow(gen_power_cpp(G_pn,n)/G_pn,1-pur) % el_pow(dgen_power_cpp(G_pn,n),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(G_pn,1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(G_dn,pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_estM_sine_cpp(mat x,mat pur,mat attri,mat price,vec beta,vec theta,mat tau_mat,int nalt,int nbeta,int ntheta){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 /(x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat ggg = el_pow(gen_sine_cpp(mylogit_cpp(G))/mylogit_cpp(G),1-pur) % el_pow(dgen_sine_cpp(mylogit_cpp(G)),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(mylogit_cpp(G),1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(mylogitpdf_cpp(G),pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_estM_sine_cpp_norm(mat x,mat pur,mat attri,mat price,vec beta,vec theta,mat tau_mat,int nalt,int nbeta,int ntheta){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 /(x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat G_dn = mat(nobs,nalt);
  mat G_pn = mat(nobs,nalt);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G_dn(i,j) = R::dnorm(G(i,j),0.0,1.0,false);
      G_pn(i,j) = R::pnorm(G(i,j),0.0,1.0,true,false);
    }
  }

  mat ggg = el_pow(gen_sine_cpp(G_pn)/G_pn,1-pur) % el_pow(dgen_sine_cpp(G_pn),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(G_pn,1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(G_dn,pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_estM_indi_cpp(mat x,mat pur,mat attri,mat price,vec beta,vec theta,int nalt,int nbeta,int ntheta){


  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  //mat ret = mylogit_cpp(G);
  mat ret = row_prods_mat(el_pow(mylogit_cpp(G),1-pur)) % J % row_prods_mat(el_pow(mylogitpdf_cpp(G),pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_estM_indi_cpp_norm(mat x,mat pur,mat attri,mat price,vec beta,vec theta,int nalt,int nbeta,int ntheta){


  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat G_dn = mat(nobs,nalt);
  mat G_pn = mat(nobs,nalt);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G_dn(i,j) = R::dnorm(G(i,j),0.0,1.0,false);
      G_pn(i,j) = R::pnorm(G(i,j),0.0,1.0,true,false);
    }
  }

  mat ret = row_prods_mat(el_pow(G_pn,1-pur)) % J % row_prods_mat(el_pow(G_dn,pur));
  return ret;

}



// [[Rcpp::export]]
mat likelihood_setM_power_cpp(mat x,mat pur,vec z,mat attri,mat price,vec beta,mat tau_mat,int nalt,int nbeta,double n){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  //gamma
  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));
  mat psi = exp(attri * beta_sub);

  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);

  //gamma
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  //gamma
  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z+1)));

  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  //gamma
  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;

  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z+1)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat ggg = el_pow(gen_power_cpp(mylogit_cpp(G),n)/mylogit_cpp(G),1-pur) % el_pow(dgen_power_cpp(mylogit_cpp(G),n),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(mylogit_cpp(G),1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(mylogitpdf_cpp(G),pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_setM_power_cpp_norm(mat x,mat pur,vec z,mat attri,mat price,vec beta,mat tau_mat,int nalt,int nbeta,double n){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;

  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat G_dn = mat(nobs,nalt);
  mat G_pn = mat(nobs,nalt);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G_dn(i,j) = R::dnorm(G(i,j),0.0,1.0,false);
      G_pn(i,j) = R::pnorm(G(i,j),0.0,1.0,true,false);
    }
  }

  mat ggg = el_pow(gen_power_cpp(G_pn,n)/G_pn,1-pur) % el_pow(dgen_power_cpp(G_pn,n),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(G_pn,1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(G_dn,pur));
  return ret;

}


// [[Rcpp::export]]
mat likelihood_setM_sine_cpp(mat x,mat pur,vec z,mat attri,mat price,vec beta,mat tau_mat,int nalt,int nbeta){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat ggg = el_pow(gen_sine_cpp(mylogit_cpp(G))/mylogit_cpp(G),1-pur) % el_pow(dgen_sine_cpp(mylogit_cpp(G)),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(mylogit_cpp(G),1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(mylogitpdf_cpp(G),pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_setM_sine_cpp_norm(mat x,mat pur,vec z,mat attri,mat price,vec beta,mat tau_mat,int nalt,int nbeta){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat G_dn = mat(nobs,nalt);
  mat G_pn = mat(nobs,nalt);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G_dn(i,j) = R::dnorm(G(i,j),0.0,1.0,false);
      G_pn(i,j) = R::pnorm(G(i,j),0.0,1.0,true,false);
    }
  }

  mat ggg = el_pow(gen_sine_cpp(G_pn)/G_pn,1-pur) % el_pow(dgen_sine_cpp(G_pn),pur);
  mat AAA = ggg * tau_mat * ggg.t();

  mat ret = row_prods_mat(el_pow(G_pn,1-pur)) % (1+AAA.diag()) % J % row_prods_mat(el_pow(G_dn,pur));
  return ret;


}




// [[Rcpp::export]]
mat likelihood_setM_indi_cpp(mat x,mat pur,vec z,mat attri,mat price,vec beta,int nalt,int nbeta){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z+1)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  //mat ret = mylogit_cpp(G);
  mat ret = row_prods_mat(el_pow(mylogit_cpp(G),1-pur)) % J % row_prods_mat(el_pow(mylogitpdf_cpp(G),pur));
  return ret;

}

// [[Rcpp::export]]
mat likelihood_setM_indi_cpp_norm(mat x,mat pur,vec z,mat attri,mat price,vec beta,int nalt,int nbeta){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  mat G_dn = mat(nobs,nalt);
  mat G_pn = mat(nobs,nalt);

  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G_dn(i,j) = R::dnorm(G(i,j),0.0,1.0,false);
      G_pn(i,j) = R::pnorm(G(i,j),0.0,1.0,true,false);
    }
  }

  mat ret = row_prods_mat(el_pow(G_pn,1-pur)) % J % row_prods_mat(el_pow(G_dn,pur));
  return ret;



}



// [[Rcpp::export]]
double sd_r_cpp_call(const Rcpp::NumericVector& x){

  // Obtain environment containing function
  Rcpp::Environment base("package:stats");

  // Make function callable from C++
  Rcpp::Function sd_r = base["sd"];

  // Call the function and receive its list output
  Rcpp::NumericVector res = sd_r(Rcpp::_["x"] = x,
                                 Rcpp::_["na.rm"]  = true); // example of additional param

  // Return test object in list structure
  return res[0];
}

// [[Rcpp::export]]
double pmvnorm_r_cpp_call(vec lower,vec upper,vec mean,mat sigma){
  //double pmvnorm_r_cpp_call(const Rcpp::NumericVector& lower,const Rcpp::NumericVector& upper,const Rcpp::NumericVector& mean,const Rcpp::NumericMatrix& sigma){

  Rcpp::Environment base("package:mvtnorm");
  Rcpp::Function pmvnorm_r = base["pmvnorm"];
  Rcpp::NumericVector res = pmvnorm_r(Rcpp::_["lower"] = lower,
                                      Rcpp::_["upper"]  = upper,
                                      Rcpp::_["mean"]  = mean,
                                      Rcpp::_["sigma"]  = sigma
                                      ); // example of additional param

  return res[0];

}

// [[Rcpp::export]]
double dmvnorm_r_cpp_call(const Rcpp::NumericVector& x,const Rcpp::NumericVector& mean,const Rcpp::NumericMatrix& sigma){

  //double dmvnorm_r_cpp_call(vec x,vec mean,mat sigma){

  Rcpp::Environment base("package:mvtnorm");
  Rcpp::Function dmvnorm_r = base["dmvnorm"];
  Rcpp::NumericVector res = dmvnorm_r(Rcpp::_["x"] = x,
                                      Rcpp::_["mean"]  = mean,
                                      Rcpp::_["sigma"]  = sigma,
                                      Rcpp::_["log"]  = false
                                      ); // example of additional param
  return res[0];

}


// [[Rcpp::export]]
vec likelihood_estM_gau_cpp(mat x,mat pur,mat attri,mat price,vec beta,vec theta,mat temp_cov,int nalt,int nbeta,int ntheta,List pur_l,List npur_l){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));
  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1/ (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  vec ret = zeros<vec>(nobs);
  vec temp_mean = zeros<vec>(nalt);
  vec temp_above = ones<vec>(nalt);

  Rcpp::Environment base("package:bayesm");
  Rcpp::Function ghkvec_r = base["ghkvec"];
  Rcpp::Function lndMvn_r = base["lndMvn"];
  
  for (int i = 0; i < nobs; i++) {
      if(sum_na_vec(pur.row(i).t())==nalt){
         ret(i) = exp(as<double>(wrap(lndMvn_r(G.row(i).t(),temp_mean,solve(chol(temp_cov), eye<mat>(nalt, nalt)))))) * J(i);
      }else if(sum_na_vec(pur.row(i).t())==0){
         ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov).t(),G.row(i).t(),temp_above,1000)));
      }else{
        int temp_pur = sum_na_vec(pur.row(i).t());
        int temp_npur = nalt - temp_pur;
        vec temp_mean_pur = zeros<vec>(temp_pur);
        vec temp_above_sub = ones<vec>(temp_npur);
        
        uvec temp_pur_l = pur_l[i];
        uvec temp_npur_l = npur_l[i];
        
        mat sigma11 = temp_cov.submat(temp_npur_l,temp_npur_l);
        mat sigma12 = temp_cov.submat(temp_npur_l,temp_pur_l);
        mat sigma21 = temp_cov.submat(temp_pur_l,temp_npur_l);
        mat sigma22 = temp_cov.submat(temp_pur_l,temp_pur_l);
        mat isigma22 = solve(sigma22, eye<mat>(temp_pur,temp_pur));
        
        uvec temp_i(1);
        temp_i(0) = i;
        vec temp_mean_sub = sigma12 * isigma22 * G(temp_i,temp_pur_l).t();
        mat temp_cov_sub = sigma11 - sigma12 * isigma22 * sigma21;
        
        //ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i);
        ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i)*exp(as<double>(wrap(lndMvn_r(G(temp_i,temp_pur_l).t(),temp_mean_pur,solve(chol(sigma22), eye<mat>(temp_pur, temp_pur))))));

      }

  }
  
  return ret;

}



// [[Rcpp::export]]
vec likelihood_estM_gau_cpp_logit(mat x,mat pur,mat attri,mat price,vec beta,vec theta,mat temp_cov,int nalt,int nbeta,int ntheta,List pur_l,List npur_l){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));

  mat psi = exp(attri * beta_sub);
  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);
  mat gamma_mat = kron(onev,trans(gamma));

  double M = exp(theta(ntheta-1));
  vec z = M - row_sums(price % x);

  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z+1)));
  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1/ (x + 1) % pur;
  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);


  mat FG = mylogit_cpp(G);
  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G(i,j) = R::qnorm(FG(i,j),0.0,1.0,true,false);

    }
  }

  vec ret = zeros<vec>(nobs);
  vec temp_mean = zeros<vec>(nalt);
  vec temp_above = ones<vec>(nalt);

  Rcpp::Environment base("package:bayesm");
  Rcpp::Function ghkvec_r = base["ghkvec"];
  Rcpp::Function lndMvn_r = base["lndMvn"];

  for (int i = 0; i < nobs; i++) {
      if(sum_na_vec(pur.row(i).t())==nalt){
         ret(i) = exp(as<double>(wrap(lndMvn_r(G.row(i).t(),temp_mean,solve(chol(temp_cov), eye<mat>(nalt, nalt)))))) * J(i);
      }else if(sum_na_vec(pur.row(i).t())==0){
         ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov).t(),G.row(i).t(),temp_above,1000)));
      }else{
        int temp_pur = sum_na_vec(pur.row(i).t());
        int temp_npur = nalt - temp_pur;
        vec temp_mean_pur = zeros<vec>(temp_pur);
        vec temp_above_sub = ones<vec>(temp_npur);

        uvec temp_pur_l = pur_l[i];
        uvec temp_npur_l = npur_l[i];

        mat sigma11 = temp_cov.submat(temp_npur_l,temp_npur_l);
        mat sigma12 = temp_cov.submat(temp_npur_l,temp_pur_l);
        mat sigma21 = temp_cov.submat(temp_pur_l,temp_npur_l);
        mat sigma22 = temp_cov.submat(temp_pur_l,temp_pur_l);
        mat isigma22 = solve(sigma22, eye<mat>(temp_pur,temp_pur));

        uvec temp_i(1);
        temp_i(0) = i;
        vec temp_mean_sub = sigma12 * isigma22 * G(temp_i,temp_pur_l).t();
        mat temp_cov_sub = sigma11 - sigma12 * isigma22 * sigma21;

        //ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i);
        ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i)*exp(as<double>(wrap(lndMvn_r(G(temp_i,temp_pur_l).t(),temp_mean_pur,solve(chol(sigma22), eye<mat>(temp_pur, temp_pur))))));

      }

  }

  return ret;

}



// [[Rcpp::export]]
vec likelihood_setM_gau_cpp(mat x,mat pur,vec z,mat attri,mat price,vec beta,mat temp_cov,int nalt,int nbeta,List pur_l,List npur_l){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  //gamma
  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));
  mat psi = exp(attri * beta_sub);


  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);

  //gamma
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  //gamma
  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  //gamma
  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;

  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);


  vec ret = zeros<vec>(nobs);
  vec temp_mean = zeros<vec>(nalt);
  vec temp_above = ones<vec>(nalt);

  Rcpp::Environment base("package:bayesm");
  Rcpp::Function ghkvec_r = base["ghkvec"];
  Rcpp::Function lndMvn_r = base["lndMvn"];

  for (int i = 0; i < nobs; i++) {
      if(sum_na_vec(pur.row(i).t())==nalt){
         ret(i) = exp(as<double>(wrap(lndMvn_r(G.row(i).t(),temp_mean,solve(chol(temp_cov), eye<mat>(nalt, nalt)))))) * J(i);
      }else if(sum_na_vec(pur.row(i).t())==0){
         ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov).t(),G.row(i).t(),temp_above,1000)));
      }else{
        int temp_pur = sum_na_vec(pur.row(i).t());
        int temp_npur = nalt - temp_pur;
        vec temp_mean_pur = zeros<vec>(temp_pur);
        vec temp_above_sub = ones<vec>(temp_npur);

        uvec temp_pur_l = pur_l[i];
        uvec temp_npur_l = npur_l[i];

        mat sigma11 = temp_cov.submat(temp_npur_l,temp_npur_l);
        mat sigma12 = temp_cov.submat(temp_npur_l,temp_pur_l);
        mat sigma21 = temp_cov.submat(temp_pur_l,temp_npur_l);
        mat sigma22 = temp_cov.submat(temp_pur_l,temp_pur_l);
        mat isigma22 = solve(sigma22, eye<mat>(temp_pur,temp_pur));

        uvec temp_i(1);
        temp_i(0) = i;
        vec temp_mean_sub = sigma12 * isigma22 * G(temp_i,temp_pur_l).t();
        mat temp_cov_sub = sigma11 - sigma12 * isigma22 * sigma21;

        //ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i);
        ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i)*exp(as<double>(wrap(lndMvn_r(G(temp_i,temp_pur_l).t(),temp_mean_pur,solve(chol(sigma22), eye<mat>(temp_pur, temp_pur))))));

      }

  }

  return ret;

}

// [[Rcpp::export]]
vec likelihood_setM_gau_cpp_logit(mat x,mat pur,vec z,mat attri,mat price,vec beta,mat temp_cov,int nalt,int nbeta,List pur_l,List npur_l){

  int nobs = x.n_rows;
  mat onev = ones<mat>(nobs,1);
  mat onev2 = ones<mat>(nalt,1);

  //gamma
  vec beta_sub = beta(span(0,nalt-1));
  vec gamma = exp(beta(span(nalt,nbeta-1)));
  mat psi = exp(attri * beta_sub);


  //mat psi = exp(attri * beta);

  mat psi_mat = reshape(psi,nalt,nobs);

  //gamma
  mat gamma_mat = kron(onev,trans(gamma));

  //double M = exp(theta(ntheta-1));
  //vec z = M - row_sums(price % x);

  //gamma
  mat G = - log(psi_mat.t()) + log(gamma_mat % x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  //mat G = - log(psi_mat.t()) + log(x + 1) + log(price) + kron(trans(onev2),log(1/(z)));

  //gamma
  mat J1 = gamma_mat / (gamma_mat % x + 1) % pur;
  //mat J1 = 1 / (x + 1) % pur;

  J1.elem(find(J1 == 0)).ones();
  vec prodJ1 = row_prods(J1);

  mat J2 = price % kron(trans(onev2),1/(z)) % pur;
  vec sumJ2 = row_sums(J2/J1);

  mat J = prodJ1 % (1 + sumJ2);

  //
  mat FG = mylogit_cpp(G);
  for(int i=0;i<nobs;i++){
    for(int j=0;j<nalt;j++){
      G(i,j) = R::qnorm(FG(i,j),0.0,1.0,true,false);

    }
  }
  //

  vec ret = zeros<vec>(nobs);
  vec temp_mean = zeros<vec>(nalt);
  vec temp_above = ones<vec>(nalt);

  Rcpp::Environment base("package:bayesm");
  Rcpp::Function ghkvec_r = base["ghkvec"];
  Rcpp::Function lndMvn_r = base["lndMvn"];

  for (int i = 0; i < nobs; i++) {
      if(sum_na_vec(pur.row(i).t())==nalt){
         ret(i) = exp(as<double>(wrap(lndMvn_r(G.row(i).t(),temp_mean,solve(chol(temp_cov), eye<mat>(nalt, nalt)))))) * J(i);
      }else if(sum_na_vec(pur.row(i).t())==0){
         ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov).t(),G.row(i).t(),temp_above,1000)));
      }else{
        int temp_pur = sum_na_vec(pur.row(i).t());
        int temp_npur = nalt - temp_pur;
        vec temp_mean_pur = zeros<vec>(temp_pur);
        vec temp_above_sub = ones<vec>(temp_npur);

        uvec temp_pur_l = pur_l[i];
        uvec temp_npur_l = npur_l[i];

        mat sigma11 = temp_cov.submat(temp_npur_l,temp_npur_l);
        mat sigma12 = temp_cov.submat(temp_npur_l,temp_pur_l);
        mat sigma21 = temp_cov.submat(temp_pur_l,temp_npur_l);
        mat sigma22 = temp_cov.submat(temp_pur_l,temp_pur_l);
        mat isigma22 = solve(sigma22, eye<mat>(temp_pur,temp_pur));

        uvec temp_i(1);
        temp_i(0) = i;
        vec temp_mean_sub = sigma12 * isigma22 * G(temp_i,temp_pur_l).t();
        mat temp_cov_sub = sigma11 - sigma12 * isigma22 * sigma21;

        //ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i);
        ret(i) = as<double>(wrap(ghkvec_r(chol(temp_cov_sub).t(),G(temp_i,temp_npur_l).t()-temp_mean_sub,temp_above_sub,1000)))*J(i)*exp(as<double>(wrap(lndMvn_r(G(temp_i,temp_pur_l).t(),temp_mean_pur,solve(chol(sigma22), eye<mat>(temp_pur, temp_pur))))));

      }

  }

  return ret;

}



