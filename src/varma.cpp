#include <vector>
#include "Rcpp.h"

using namespace std;
using namespace Rcpp;


class CMatrix {
public:
   std::vector< std::vector<double> > elements;
   
   CMatrix(double ele, int n_row, int n_col);
   CMatrix(std::vector< std::vector<double> > eles);
   CMatrix();
   
   int nrow();
   int ncol();
   
   double * operator()(int p, int q);
   std::vector<double> operator() (int p, bool is_row);
   
   void append(std::vector<double>);
   void append(CMatrix);
   
   double element_sum();
private:
   
};

int CMatrix::nrow(){
   return (int)elements.size();
}

int CMatrix::ncol(){
   int size;
   if (elements.size() > 0) {
      size = (int) elements[0].size();
   } else {
      size = 0;
   }
   
   return size;
}


CMatrix::CMatrix(std::vector< std::vector<double> > eles) {
   elements = eles;
}

CMatrix::CMatrix(double ele, int n_row, int n_col) {
   for (int i = 0; i < n_row; i++){
      std::vector<double> temp_col(n_col);
      fill(temp_col.begin(), temp_col.end(),ele);
      elements.push_back(temp_col);
   }
}

CMatrix::CMatrix(){
}

double CMatrix::element_sum() {
   double sum = 0;
   for (int i = 0; i < nrow(); i++){
      for (int j = 0; j < ncol(); j++) {
         sum += elements[i][j];
      }
   }
   return sum;
}

void CMatrix::append(CMatrix M){
   if (ncol() == 0 || ncol() == M.ncol()){
      for (int i = 0; i< M.nrow();i++) {
         elements.push_back(M.elements[i]);
      }
   } else {
      std::cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
}


void CMatrix::append(std::vector<double> row){
   if (ncol() == 0 || ncol() == row.size()){
      elements.push_back(row);
   } else {
      std::cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
}


double * CMatrix::operator()(int p, int q){
   return & elements[p-1][q-1];
}

std::vector<double> CMatrix::operator()(int p, bool is_row = true) {
   std::vector<double> ret;
   if (is_row) {
      ret = elements[p-1];
   } else {
      for (int i = 0; i < elements.size(); i++){
         ret.push_back(elements[i][p-1]);
      }
   }
   
   return ret;
}


CMatrix ToCMatrix (Rcpp::NumericMatrix RMat){
   CMatrix C;
   int nrow = RMat.nrow();
   int ncol = RMat.ncol();
   for (int i = 0; i < nrow; i++){
      std::vector<double> row;
      for (int j = 0; j < ncol; j++) {
         row.push_back(RMat(i,j));
      }
      C.append(row);
   }
   
   return C;
}

//Rcpp::NumericMatrix ToMatrix (CMatrix M){
//   Rcpp::NumericMatrix C(M.nrow(), M.ncol());
//   for (int i = 0; i < M.nrow(); i++) {
//      for (int j =0; j < M.ncol(); j++) {
//         C(i,j) = *M(i+1,j+1);
//      }
//   }
//   return C;
//}

CMatrix rbind(CMatrix A, CMatrix  B) {
   CMatrix C;
   if (A.ncol()==B.ncol()){
      C = A;
      for (int i = 0; i< B.nrow();i++) {
         C.elements.push_back(B.elements[i]);
      }
   } else {
      std::cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
   
   return C;
}

CMatrix as_matrix (std::vector<double> x, bool is_as_row = true){
   CMatrix C;
   if(is_as_row){
      C.elements.push_back(x);
   } else {
      for (int i = 0; i < x.size(); i++){
         std::vector<double> t;
         t.push_back(x[i]);
         C.elements.push_back(t);
      }
   }
   
   return C;
}


CMatrix prod(CMatrix A, CMatrix B) {
   CMatrix C(0, A.nrow(),B.ncol());
   if (A.ncol() == B.nrow()) {
      for (int i = 0; i < A.nrow(); i++) {
         for (int j = 0; j < B.ncol(); j++) {
            double temp = 0;
            
            for (int k = 0 ; k < B.nrow(); k++) {
               temp += (*A(i+1,k+1)) * (*B(k+1,j+1));
            }
            
            C.elements[i][j]= temp;
         }
      }
   } else {
      std::cerr<<"CMatrix Prod Error:  Incompatible matrices.\n";
   }
   
   return C;
}



class Varma {
public:
   CMatrix Obs;
   CMatrix Residuals;
   int k;  // # Param
   int nT; // # Obs
   
   std::vector<double> Ph0; // Constants
   CMatrix PH; // AR Coeff
   CMatrix TH; // MA Coeff
   
   int p;  // AR order p
   int q;  // MA order q
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() ==0
   
   Varma(CMatrix &, CMatrix &, std::vector<double> &, int, int, bool);
   
private:
   int checkMaskFormat(CMatrix &);
   void fillParamFixed(CMatrix &, std::vector<double>, bool);
   void compResiduals();
};

int Varma::checkMaskFormat(CMatrix & Mask) {
   int sum = 0;
   for (int i = 1; i <= Mask.nrow(); i++){
      for (int j = 1; j<= Mask.ncol(); j++) {
         if (*Mask(i,j)==1){
            sum += 1;
         } else if (*Mask(i,j) == 0) {
            // Skip zero elements
         } else {
            std::cerr<<"Invalid Mask CMatrix:  CMatrix contains elements others than 0 or 1 ("
            <<i<<", "<<j<<") = "<<*Mask(i,j)<<"\n";
            break;
         }
      }
   }
   return sum;
}


void Varma::fillParamFixed(CMatrix & Mask, std::vector<double> ParamFixed, bool isMeanIncluded) {
   CMatrix Beta;
   int kp = k * p;
   int kq = k * q;
   
   int i_start = 0;
   
   std::vector<double> QParamFixed(ParamFixed.size());
   reverse_copy(ParamFixed.begin(), ParamFixed.end(), QParamFixed.begin());
   
   if (hasMean){
      // Initiate Ph0 to be zero
      Ph0.resize(k);
      fill(Ph0.begin(),Ph0.end(),0);
      for (int i = 1; i <= k; i++) {
         if (*Mask(1,i) == 1) {
            Ph0.at(i-1) = QParamFixed.back();
            QParamFixed.pop_back();
         }
      }
      i_start = 1;
   }
   
   if (p > 0){
      PH = CMatrix(0, kp, k);
      for (int j = 1; j <= k; j++) {
         for (int i = 1; i<= kp; i++) {
            if (*Mask(i_start+i, j) == 1) {
               PH.elements[i-1][j-1] = QParamFixed.back();
               QParamFixed.pop_back();
            }
         }
      }
      i_start += kp;
   }
   
   if (q > 0){
      TH = CMatrix(0, kq, k);
      for (int j = 1; j<= k; j++) {
         for (int i = 1; i <= kq; i++) {
            if (*Mask(i_start+i,j) == 1) {
//               *TH(i_start+i, j) = QParamFixed.back();
               TH.elements[i-1][j-1] = QParamFixed.back();
               QParamFixed.pop_back();
            }
         }
      }
   } else if (QParamFixed.size() != 0){
      std::cerr<<"Init with parameters error:  Too many parameters\n";
   }

}


void Varma::compResiduals(){
   // Const Row
   std::vector<double> Res_Const;
   for (int i = 1; i <= k; i++){
      Res_Const.push_back(*Obs(1,i)-Ph0[i-1]);
   }
   
   Residuals.append(Res_Const);

   // ------------------------------------------------------------------------
   // Step 1:
   // Calculate std::max(p,q) residuals, necessary for subsequent recursions
   if (std::max(p,q)>1){
      for (int t = 2; t<= std::max(p,q); t++){
         std::vector<double> Res_Row;
         for (int i = 1; i <= k; i++){
            Res_Row.push_back(*Obs(t,i)-Ph0[i-1]);
         }
         
         for (int j = 1; j <= p; j++){
            if (t-j > 0){
               CMatrix PH_Slice;
               int jdx = (j-1)*k;
               for (int r = 1; r<= k; r++){
                  PH_Slice.elements.push_back(PH(jdx+r));
               }
               
               // Downcasting matrix to row std::vector
               std::vector<double> Estimate_AR = prod(as_matrix(Obs(t-j)),PH_Slice)(1);
               
               for (int r = 0; r < k; r++){
                  Res_Row[r] = Res_Row[r]-Estimate_AR[r];
               }
               
            }
         }
         
         for (int j = 1; j <= q; j++){
            if (t-j > 0){
               int jdx = (j-1)*k;
               CMatrix TH_Slice;
               for (int r = 1; r<=k; r++){
                  TH_Slice.elements.push_back(TH(jdx+r));
               }
               
               // Downcasting matrix to row std::vector
               std::vector<double> Estimate_MA = prod(as_matrix(Residuals(t-j)),TH_Slice)(1);
               
               for (int r = 0; r < k; r++){
                  Res_Row[r] = Res_Row[r]-Estimate_MA[r];
               }
            }
         }
         Residuals.append(Res_Row);
      }
   }
   
   // ------------------------------------------------------------------------
   // Step 2:
   // Calculate Residuals from Index std::max(p,q) onwards
   CMatrix Beta;
   Beta.append(Ph0);
   Beta.append(PH);
   Beta.append(TH);
   
//   // Debugging Print Beta
//   cout<<"Beta\n";
//   for (int i = 0; i < Beta.nrow(); i++){
//      copy(Beta.elements[i].begin(), Beta.elements[i].end(), ostream_iterator<double>(cout,","));
//      cout<<"\n";
//   }
   
   double Obs_Const;
   if (hasMean){
      Obs_Const = 1;
   }
   
   int i_start = std::max(p,q)+1;
   for (int t = i_start; t<= nT; t++) {
      std::vector<double> Obs_Past;
      if (p > 0) {
         for (int j=1; j<=p; j++){
            std::vector<double> Zt_Slice = Obs(t-j);
            Obs_Past.insert(Obs_Past.end(), Zt_Slice.begin(),Zt_Slice.end());
         }
      }
      
      if (q > 0) {
         for (int j=1; j<=q; j++){
            std::vector<double> At_Slice = Residuals(t-j);
            Obs_Past.insert(Obs_Past.end(), At_Slice.begin(),At_Slice.end());
         }
      }
      
      Obs_Past.insert(Obs_Past.begin(), Obs_Const);

      std::vector<double> Estimate_ARMA = prod(as_matrix(Obs_Past), Beta)(1);

      std::vector<double> Res_Row;
      for (int r = 0; r < Estimate_ARMA.size(); r++){
         Res_Row.push_back(*Obs(t,r+1)-Estimate_ARMA[r]);
      }
      Residuals.append(Res_Row);
   }
   
   // ------------------------------------------------------------------------
   // Step 3:
   // Erase the first std::max(p,q) residuals
   Residuals.elements.erase (Residuals.elements.begin(),
                             Residuals.elements.begin() + i_start-1);
}



// Model specification:  PH * Obs = PH0 + TH * AT
// FIX1:  Boolean matrix, comparing with ParamFixed, indicating fixed parameters
Varma::Varma (CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed, int ar_p, int ma_q, bool isMeanIncluded) {
   Obs = TimeSeries;
   
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   
   p = ar_p;
   q = ma_q;
   
   // Check compactible:  Mask and Parameters
   int mask_size = checkMaskFormat(Mask);
   
   if (mask_size != ParamFixed.size()) {
      std::cerr<<"Input error:  Fixed mask and parameters provided are of different size.\n";
      std::cerr<<"   Mask contains "<<mask_size<<" 1s but "<<ParamFixed.size()<<" is provided.";
   } else {
      // Check compactible:  Mask and p, q, k
      int row_const = hasMean? 1 : 0;
      
      if (row_const+k*(p+q) != Mask.nrow()) {
         std::cerr<<"Input error:  Fixed mask has a wrong row number.\n";
         std::cerr<<"Need "<<row_const+k*(p+q)<<" rows but the mask matrix has "<<Mask.nrow()<<" rows.\n";
      } else {
         fillParamFixed(Mask, ParamFixed, hasMean);
         
         compResiduals();
      }
   }
}



RcppExport SEXP GetVarmaResiduals(SEXP _timeSeries,
                                  SEXP _mask,
                                  SEXP _paramFixed,
                                  SEXP _p,
                                  SEXP _q,
                                  SEXP _isMeanIncluded
                                  ) {
   
   Rcpp::NumericMatrix RTimeSeries(_timeSeries);
   Rcpp::NumericMatrix RMask(_mask);
   
   CMatrix TimeSeries = ToCMatrix(RTimeSeries);
   CMatrix Mask = ToCMatrix(RMask);
   
   std::vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< std::vector<double> >(_paramFixed);
   }
   
   int ar_p = Rcpp::as<int>(_p);
   int ma_q = Rcpp::as<int>(_q);
   bool isMeanIncluded = Rcpp::as<bool>(_isMeanIncluded);

   Varma varma(TimeSeries, Mask, ParamFixed, ar_p, ma_q, isMeanIncluded);
   
   
   CMatrix Residuals = varma.Residuals;
   return Rcpp::wrap(Residuals.elements);
}


