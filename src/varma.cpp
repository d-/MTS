#include <vector>
#include <iostream>

#include "Rcpp.h"

using namespace std;
using namespace Rcpp;

// MARK: CMatrix
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
   void transpose();
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
      //      cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
}

void CMatrix::append(std::vector<double> row){
   if (ncol() == 0 || ncol() == row.size()){
      elements.push_back(row);
   } else {
      //      cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
}

void CMatrix::transpose() {
   std::vector< std::vector<double> > transposed_elements;
   if (elements.size() > 0){
      for (int j = 0; j < elements[0].size(); j++) {
         std::vector<double> temp;
         for (int i = 0; i < elements.size(); i++) {
            temp.push_back(elements[i][j]);
         }
         transposed_elements.push_back(temp);
      }
   }
   
   elements = transposed_elements;
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


// MARK: CMatrix Utilities

CMatrix Cnegative(CMatrix &A) {
   CMatrix B = A;
   if (A.elements.size() > 0){
      for (int i = 0; i < A.elements.size(); i++) {
         for (int j = 0; j < A.elements[0].size(); j++) {
            B.elements[i][j] = - A.elements[i][j];
         }
      }
   }
   return (B);
}

CMatrix Ctranspose(CMatrix &A) {
   CMatrix B;
   if (A.elements.size() > 0){
      for (int j = 0; j < A.elements[0].size(); j++) {
         std::vector<double> temp;
         for (int i = 0; i < A.elements.size(); i++) {
            temp.push_back(A.elements[i][j]);
         }
         B.elements.push_back(temp);
      }
   }
   return (B);
}

CMatrix Cdiagonal (std::vector<double> &A) {
   CMatrix B (0, (int) A.size(), (int) A.size());
   for (int i = 0; i < A.size(); i++){
      B.elements[i][i] = A[i];
   }
   return (B);
}

CMatrix Cidentity (int s){
   std::vector<double> ones(s);
   fill(ones.begin(), ones.end(), 1);
   CMatrix A = Cdiagonal(ones);
   return (A);
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

CMatrix rbind(CMatrix A, CMatrix  B) {
   CMatrix C;
   if (A.ncol()==B.ncol()){
      C = A;
      for (int i = 0; i< B.nrow();i++) {
         C.elements.push_back(B.elements[i]);
      }
   } else {
      //      cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
   
   return C;
}

CMatrix cbind(CMatrix A, CMatrix  B) {
   CMatrix C;
   if (A.nrow()==B.nrow()){
      C = A;
      for (int i = 0; i< A.nrow();i++) {
         C.elements[i].insert(C.elements[i].end(), B.elements[i].begin(), B.elements[i].end());
      }
   } else {
      //      cerr<<"rbind Error:  Can't cbind matrices with different row sizes.\n";
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
      //      cerr<<"CMatrix Prod Error:  Incompatible matrices.\n";
   }
   
   return C;
}

CMatrix rows(CMatrix &A, int i0, int i1) {
   CMatrix C;
   for (int i = i0; i <= i1; i++) {
      C.append(A(i));
   }
   return C;
}

CMatrix cols(CMatrix &A, int i0, int i1) {
   CMatrix C;
   for (int i = i0; i <= i1; i++) {
      C.append(A(i, false));
   }
   
   C.transpose();
   return C;
}


CMatrix matrix_prod (CMatrix &A, CMatrix &B, int p, int P)
{
   int k = A.nrow();
   CMatrix C = A;
   for (int i = 1; i <= P; i++) {
      int i_start = (i-1) * k;
      CMatrix m2 (cols(B,i_start+1, i_start+k));
      C = cbind(C,m2);
      for (int j  = 1; j <= p; j++){
         int j_start = (j-1) * k;
         CMatrix m1 (cols(A,j_start+1, j_start+k));
         C= cbind(C, prod(Cnegative(m1),m2));
      }
   }
   
   return C;
}

CMatrix matrix_prod_alt (CMatrix &A, CMatrix &B, int p, int P)
{
   int k = A.nrow();
   CMatrix C = A;
   for (int i = 1; i <= P; i++) {
      int i_start = (i-1) * k;
      CMatrix m2 (cols(B,i_start+1, i_start+k));
      C = cbind(C,m2);
      for (int j  = 1; j <= p; j++){
         int j_start = (j-1) * k;
         CMatrix m1 (cols(A,j_start+1, j_start+k));
         C= cbind(C, prod(Cnegative(m2),m1));
      }
   }
   return C;
}


// MARK: Varma

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
            //            cerr<<"Invalid Mask CMatrix:  CMatrix contains elements others than 0 or 1 ("
            //            <<i<<", "<<j<<") = "<<*Mask(i,j)<<"\n";
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
      for (int i = 1; i <= k; i++) {
         if (*Mask(1,i) == 1) {
            Ph0.at(i-1) = QParamFixed.back();
            QParamFixed.pop_back();
         }
      }
      i_start = 1;
   }
   
   if (p > 0){
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
      //      cerr<<"Init with parameters error:  Too many parameters\n";
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
   // Calculate max(p,q) residuals, necessary for subsequent recursions
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
   // Calculate Residuals from Index max(p,q) onwards
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
   // Erase the first max(p,q) residuals
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
   
   // Parametric initializations
   Ph0.resize(k);
   fill(Ph0.begin(),Ph0.end(),0);
   PH = CMatrix(0, k*p, k);
   TH = CMatrix(0, k*q, k);
   
   if (ParamFixed.size()>0){
      fillParamFixed(Mask, ParamFixed, hasMean);
   }
   
   compResiduals();
   
}

// MARK: SVarma
class SVarma {
public:
   CMatrix Obs;
   CMatrix Residuals;
   int k;  // # Param
   int nT; // # Obs
   
   std::vector<double> Ph0; // Constants
   CMatrix PH; // AR Coeff
   CMatrix TH; // MA Coeff
   CMatrix sPH; // sAR Coeff
   CMatrix sTH; // sMA Coeff
   
   bool matrix_prod_method_switch; // swi
   CMatrix Phi; //
   CMatrix Theta; //
   
   CMatrix Beta; //!!! Local property in Varma
   std::vector<int> ARlags;
   std::vector<int> MAlags;
   
   int nar;  // # Param
   int nma; // # Obs
   
   int p;  // AR order p
   int q;  // MA order q
   int P;  // SAR order P
   int Q;  // SMA order Q
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() ==0
   
   SVarma (CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed, bool isMeanIncluded, std::vector<int> Orders, std::vector<int> ARlags, std::vector<int> MAlags, CMatrix & Sresi, bool swi);
private:
   int checkMaskFormat(CMatrix & Mask);
   void fillParamFixed(CMatrix & Mask, std::vector<double> ParamFixed, bool isMeanIncluded);
   void compResiduals(CMatrix & Sresiduals);
};

int SVarma::checkMaskFormat(CMatrix & Mask) {
   int sum = 0;
   for (int i = 1; i <= Mask.nrow(); i++){
      for (int j = 1; j<= Mask.ncol(); j++) {
         if (*Mask(i,j)==1){
            sum += 1;
         } else if (*Mask(i,j) == 0) {
            // Skip zero elements
         } else {
            //            cerr<<"Invalid Mask CMatrix:  CMatrix contains elements others than 0 or 1 ("
            //            <<i<<", "<<j<<") = "<<*Mask(i,j)<<"\n";
            break;
         }
      }
   }
   return sum;
}

void SVarma::fillParamFixed(CMatrix & Mask, std::vector<double> ParamFixed, bool isMeanIncluded) {
   CMatrix Beta;
   int kp = k * p;
   int kq = k * q;
   int kP = k * P;
   int kQ = k * Q;
   
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
   
   
   if (nar > 0){
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
      
      if (P > 0){
         sPH = CMatrix(0, kP, k);
         for (int j = 1; j <= k; j++) {
            for (int i = 1; i<= kP; i++) {
               if (*Mask(i_start+i, j) == 1) {
                  sPH.elements[i-1][j-1] = QParamFixed.back();
                  QParamFixed.pop_back();
               }
            }
         }
         i_start += kP;
      }
      
   }
   
   if (nar >0){
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
         //         cerr<<"Init with parameters error:  Too many parameters\n";
      }
      
      if (Q > 0){
         sTH = CMatrix(0, kQ, k);
         for (int j = 1; j<= k; j++) {
            for (int i = 1; i <= kQ; i++) {
               if (*Mask(i_start+i,j) == 1) {
                  //               *TH(i_start+i, j) = QParamFixed.back();
                  sTH.elements[i-1][j-1] = QParamFixed.back();
                  QParamFixed.pop_back();
               }
            }
         }
      } else if (QParamFixed.size() != 0){
         //         cerr<<"Init with parameters error:  Too many parameters\n";
      }
   }
   
   // Building Beta Matrix
   Beta.append(Ph0);
   if (p>0 && p>0){
      if (matrix_prod_method_switch) {
         Phi = matrix_prod_alt(PH, sPH, p, P);
      }
      else {
         Phi = matrix_prod(PH,sPH,p,P);
      }
      Beta.append(Ctranspose(Phi));
   }
   
   if (p>0 && P==0) {
      Beta.append(Ctranspose(PH));
   }
   
   if (p==0 && P>0) {
      Beta.append(Ctranspose(sPH));
   }
   
   if (q>0 && Q>0){
      if (matrix_prod_method_switch) {
         Theta = matrix_prod_alt(TH, sTH, q, Q);
      }
      else {
         Theta = matrix_prod(TH, sTH, q, Q);
      }
      CMatrix Theta_transposed = Ctranspose(Theta);
      Beta.append(Cnegative(Theta_transposed));
   }
   
   if (q>0 && Q==0) {
      CMatrix TH_transposed = Ctranspose(TH);
      Beta.append(Cnegative(TH_transposed));
   }
   
   if (q==0 && Q>0) {
      CMatrix sTH_transposed = Ctranspose(sTH);
      Beta.append(Cnegative(sTH_transposed));
   }
}

void SVarma::compResiduals(CMatrix & Sresi) {
   Residuals = Sresi;
   
   double Obs_Const;
   if (hasMean){
      Obs_Const = 1;
   }
   
   int i_start = std::max(nar,nma)+1;
   for (int t = i_start; t<= nT; t++) {
      std::vector<double> Obs_Past;
      if (nar > 0) {
         for (int j=1; j<=nar; j++){
            int jj = ARlags[j-1];
            std::vector<double> Zt_Slice = Obs(t-jj);
            Obs_Past.insert(Obs_Past.end(), Zt_Slice.begin(),Zt_Slice.end());
         }
      }
      
      if (nma > 0) {
         for (int j=1; j<=nma; j++){
            int jj = MAlags[j-1];
            std::vector<double> At_Slice = Residuals(t-jj);
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
   
   Residuals.elements.erase (Residuals.elements.begin(),
                             Residuals.elements.begin() + i_start-1);
}

SVarma::SVarma (CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed, bool isMeanIncluded, std::vector<int> Orders, std::vector<int> ARlags, std::vector<int> MAlags, CMatrix & Sresi, bool swi)
{
   Obs = TimeSeries;
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   
   p = Orders[0];  // AR_p
   q = Orders[2];  // MA_q
   
   P = Orders[3];  // SAR_p
   Q = Orders[5];  // SMA_q
   
   nar = (int) ARlags.size();
   nma = (int) MAlags.size();
   
   int i_start = 0;
   i_start = std::max(ARlags.back(), MAlags.back())+1;
   
   Ph0.resize(k);
   fill(Ph0.begin(),Ph0.end(),0);
   PH = CMatrix(0, k*p, k);
   TH = CMatrix(0, k*q, k);
   
   if (ParamFixed.size()>0){
      fillParamFixed(Mask, ParamFixed, hasMean);
   }
   
   compResiduals(Sresi);
   
   //   // Check compactible:  Mask and Parameters
   //   int mask_size = checkMaskFormat(Mask);
   //
   //   if (ParamFixed.size() > 0 && mask_size != ParamFixed.size()) {
   //      //      cerr<<"Input error:  Fixed mask and parameters provided are of different size.\n";
   //      //      cerr<<"   Mask contains "<<mask_size<<" 1s but "<<ParamFixed.size()<<" is provided.";
   //   } else {
   //      // Check compactible:  Mask and p, q, k
   //      int row_const = hasMean? 1 : 0;
   //
   //      if (Mask.nrow()>0 && row_const+k*(nar+nma) != Mask.nrow()) {
   //         //         cerr<<"Input error:  Fixed mask has a wrong row number.\n";
   //         //         cerr<<"Need "<<row_const+k*(p+q)<<" rows but the mask matrix has "<<Mask.nrow()<<" rows.\n";
   //      } else {
   //         if (Mask.nrow()>0){
   //         }
   //      }
   //   }
   
   
}


// MARK: VMA
class VMA {
public:
   CMatrix Obs;
   CMatrix Residuals;
   int k;  // # Param
   int nT; // # Obs
   
   int q;  // MA order q
   
   std::vector<double> Ph0; // Constants -- OC: mu
   
   CMatrix Theta;
   CMatrix TH;
   
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() == 0
   
   VMA(CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed,int ma_q, bool isMeanIncluded);
private:
};

VMA::VMA(CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed,int ma_q, bool isMeanIncluded)
{
   Obs = TimeSeries;
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   q = ma_q;
   
   std::vector<double> QParamFixed(ParamFixed.size());
   reverse_copy(ParamFixed.begin(), ParamFixed.end(), QParamFixed.begin());
   
   int i_start = 0;
   
   // Initiate Ph0 to be zero
   Ph0.resize(k);
   fill(Ph0.begin(),Ph0.end(),0);
   
   if (hasMean){
      i_start = 1;
      for (int i = 1; i <= k; i++) {
         if (*Mask(1,i) == 1) {
            Ph0.at(i-1) = QParamFixed.back();
            QParamFixed.pop_back();
         }
      }
      
      // Remove the mean
      for (int j = 1; j <= k; j++) {
         for (int i = 1; i <= nT; i++){
            Obs.elements[i-1][j-1] = Obs.elements[i-1][j-1]-Ph0[j-1];
         }
      }
   }
   
   int kq  =  k*q;
   Theta = CMatrix(0, kq, k);
   for (int j = 1; j <= k; j++) {
      for (int i = 1; i<= kq; i++) {
         if (*Mask(i_start+i, j) == 1) {
            Theta.elements[i-1][j-1] = QParamFixed.back();
            QParamFixed.pop_back();
         }
      }
   }
   i_start += kq;
   
   // Theta = rbind[theta_1',theta_2', ..., theta_q']
   // Checking the invertibility of t(Theta)
   TH = Ctranspose(Theta);
   
   if (q > 1) {
      CMatrix ones = Cidentity((q-1)*k);
      CMatrix zeros(0, (q-1)*k, k);
      CMatrix tmp = cbind(ones, zeros);
      TH.append(tmp);
   }
}


class VMADemean {
public:
   CMatrix Obs;
   CMatrix Residuals;
   int k;  // # Param
   int nT; // # Obs
   
   int q;  // MA order q
   
   std::vector<double> Ph0; // Constants -- OC: mu
   
   CMatrix Theta;
   CMatrix TH;
   
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() == 0
   
   VMADemean(CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed,int ma_q, bool isMeanIncluded);
private:
};

VMADemean::VMADemean(CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed,int ma_q, bool isMeanIncluded)
{
   Obs = TimeSeries;
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   q = ma_q;
   
   std::vector<double> QParamFixed(ParamFixed.size());
   reverse_copy(ParamFixed.begin(), ParamFixed.end(), QParamFixed.begin());
   
   int i_start = 0;
   
   if (hasMean){
      i_start = 1;
      // Initiate Ph0 to be zero
      Ph0.resize(k);
      fill(Ph0.begin(),Ph0.end(),0);
      for (int i = 1; i <= k; i++) {
         if (*Mask(1,i) == 1) {
            Ph0.at(i-1) = QParamFixed.back();
            QParamFixed.pop_back();
         }
      }
      
      // Remove the mean
      for (int j = 1; j <= k; j++) {
         for (int i = 1; i <= nT; i++){
            Obs.elements[i-1][j-1] = Obs.elements[i-1][j-1]-Ph0[j-1];
         }
      }
   }
}

//// MARK: C++ Tests
//// DISABLE BEFORE RCPP COMPILATION
//// Simple Testing Hook
//#include <iostream>
//int main()
//{
//   //   Varma::Varma (CMatrix & TimeSeries, CMatrix & Mask, std::vector<double> & ParamFixed, int ar_p, int ma_q, bool isMeanIncluded)
//   vector<vector<double>> v_ts = {{-0.68025910, -1.75375261},
//      {-0.70460252, -3.03033729},
//      {-0.42229690, -3.36407270},
//      { 0.09358394, -4.73390906},
//      {-1.06262633, -4.47887795},
//      {-2.58585688, -3.27856854},
//      {-0.45325621, -1.19216546},
//      { 1.77141156,  1.35579414},
//      { 1.43689844,  2.36952145},
//      { 1.04269587,  1.64962734},
//      { 1.75173443, -0.18663963},
//      { 1.56053934, -3.14823630},
//      { 1.32138422, -4.17667064},
//      { 0.05338261, -4.82166510},
//      {-1.43369349, -3.82470517},
//      {-0.59601647, -2.30285534},
//      {-0.95794451, -3.88044996},
//      {-2.84679183, -4.45443046},
//      {-0.93553139, -3.90337400},
//      {-0.23827492, -3.26142641},
//      {-0.94046561, -3.13019564},
//      {-1.65130002, -1.37351954},
//      {-0.90747956, -2.47878531},
//      {-1.50023182, -2.93223254},
//      {-1.61499384, -3.21890847},
//      {-2.22842297, -3.73085572},
//      {-1.69707802, -3.74197848},
//      {-0.10245939, -2.87965686},
//      {-0.58907041, -2.34942438},
//      {-0.68389708, -2.01944429},
//      {-1.05598147,  0.29354952},
//      { 0.76799345,  0.26712775},
//      { 1.26318299, -1.03430634},
//      {-1.22634540, -1.13142289},
//      {-0.64085912, -0.60011329},
//      {-0.69292349, -0.70204712},
//      {-0.11714641, -0.08132956},
//      { 1.57729783,  1.70695745},
//      { 3.44119364,  1.20067348},
//      { 2.66655378, -2.27100279},
//      {-0.24805538, -5.14870120},
//      {-2.29621527, -5.42722889},
//      {-3.14357079, -3.60744809},
//      {-3.11935124, -0.80036247},
//      {-0.26093559,  3.11525904},
//      { 0.99131232,  4.81572486},
//      { 0.96802945,  4.62159388},
//      { 3.84489899,  4.40476317},
//      { 3.99733936,  2.17658044},
//      { 3.83053354, -1.31124235},
//      { 1.19575846, -2.79092074},
//      { 0.30265243, -2.91735623},
//      {-0.32930897, -2.19985726},
//      {-1.44053764, -3.90528105},
//      {-1.01967940, -6.11864947},
//      {-2.47135405, -5.87573412},
//      {-2.42696414, -3.77640441},
//      {-0.43006170, -2.90777023},
//      {-1.21733319, -2.69930121},
//      {-1.36257994, -0.54042603},
//      {-0.40483841,  1.48124274},
//      {-1.95976214,  3.52922861},
//      {-1.22304837,  7.40997150},
//      { 3.67603232, 10.96243649},
//      { 4.81774261, 11.38422439},
//      { 4.78828936, 10.24558976},
//      { 4.37502134, 10.02160477},
//      { 3.63546452, 10.28458210},
//      { 3.49982540,  9.87308236},
//      { 3.64367503,  9.08972054},
//      { 5.29157531,  8.28366292},
//      { 4.70225241,  5.58957257},
//      { 1.77210909,  2.80942612},
//      {-0.02240048,  2.67997069},
//      { 0.61835493,  1.60075577},
//      { 2.08967422,  0.48781193},
//      { 2.39928854, -0.15878012},
//      { 1.87679777, -0.80318810},
//      {-0.92807732, -3.95727100},
//      {-3.49453768, -5.68721124},
//      {-1.88224143, -4.39916762},
//      { 0.22583413, -3.58061144},
//      {-1.53914061, -3.63459269},
//      {-2.52644586, -3.24150277},
//      {-1.37236890, -2.49412715},
//      {-1.48557174, -0.88794990},
//      { 0.60119516,  0.78088283},
//      { 0.58163580, -0.37615165},
//      {-1.42839359,  0.19376096},
//      {-1.60915525,  2.54253265},
//      { 0.47475990,  5.13307899},
//      { 2.19517619,  6.15565934},
//      { 2.08476873,  5.78298467},
//      { 2.39708128,  5.23715461},
//      { 2.03865750,  2.40170755},
//      { 0.93942648, -0.94240319},
//      {-0.71285290, -1.71405182},
//      {-0.43347394, -2.43485146},
//      {-1.14504194, -3.86150678},
//      {-0.39159677, -2.32558354},
//      {-1.18325693,  0.34365839},
//      {-0.80675953,  1.74171090},
//      {-2.57429273,  3.63342451},
//      { 0.28620406,  6.25357840},
//      { 2.41352566,  6.72105974},
//      { 1.61483506,  6.53026073},
//      { 2.34896730,  6.16693173},
//      { 2.83533172,  4.38092947},
//      { 0.24461950,  3.28455944},
//      { 0.60312944,  3.35566687},
//      { 2.50133092,  2.36280439},
//      { 2.90280013,  0.93207658},
//      { 1.63865628, -0.81556168},
//      { 0.54470369, -1.03822774},
//      { 0.28909755,  0.15750390},
//      {-0.28319372,  0.51822590},
//      {-0.74126070, -0.15266919},
//      {-0.62452369, -1.43490400},
//      {-1.66175049, -0.96383439},
//      {-0.18552548,  2.31619099},
//      {-0.48912391,  2.84341817},
//      { 1.73354272,  0.97441219},
//      { 2.22026813, -2.54927177},
//      {-0.10584969, -5.58156599},
//      {-0.16180690, -7.00978865},
//      {-1.06256961, 10.12924788},
//      {-2.99490422, 11.24724314},
//      {-2.37553367, 11.96271086},
//      {-2.81602366, 12.18681783},
//      {-5.06941404, 11.02234118},
//      {-4.81285225, -9.12290615},
//      {-3.44895033, -6.56074913},
//      {-1.92871349, -3.25744003},
//      {-2.53604676, -1.81695309},
//      {-2.08160501,  1.21903070},
//      {-0.18340175,  3.77558305},
//      { 0.95423515,  4.89582324},
//      { 0.76206986,  3.72047504},
//      { 2.19031583,  1.56360928},
//      { 4.01094876,  1.69344396},
//      { 3.00494360,  0.39443698},
//      { 0.88360416, -1.44694533},
//      {-0.02336771, -2.94541803},
//      {-1.60900832, -3.61165905},
//      {-0.23458038, -3.78317660},
//      {-1.52604787, -5.17740496},
//      {-1.41426019, -4.01178882},
//      {-1.57751183, -2.96659282},
//      {-2.77902214, -2.58491762},
//      {-1.18092764, -2.41798573},
//      {-0.08484436, -2.66020165},
//      {-1.87051270, -2.12853789},
//      {-2.52979871, -1.06075363},
//      {-1.07359316, -0.06177132},
//      {-0.52942827,  0.83116242},
//      {-2.04541987,  1.96613877},
//      {-1.18043307,  4.92320698},
//      { 3.79722223,  6.71049257},
//      { 4.15093264,  4.10548404},
//      { 1.08206919,  3.49848972},
//      { 1.75148887,  4.39174050},
//      { 2.11206019,  4.02047963},
//      { 2.05900856,  1.61991040},
//      { 1.48894743,  0.41242995},
//      { 0.49915364,  1.66718411},
//      { 1.21332411,  2.99836097},
//      { 1.62015804,  3.40421605},
//      { 0.75460282,  3.16190635},
//      { 1.51026651,  3.74227411},
//      { 1.80125854,  3.58285387},
//      {-0.36416386,  3.52550979},
//      {-0.09568142,  4.04055479},
//      { 2.11225027,  5.50511978},
//      { 2.39015013,  5.23623392},
//      { 1.66012525,  2.90662043},
//      { 3.29299078,  0.67889644},
//      { 2.76338105, -1.12060006},
//      { 1.54113222, -1.93147629},
//      {-0.39132530, -3.00409535},
//      {-1.14838617, -3.58381921},
//      {-0.61454167, -2.29497468},
//      {-0.06210385, -2.20686155},
//      {-0.73570233, -2.58134869},
//      {-2.54112257, -3.87616492},
//      {-3.10442131, -4.43529887},
//      {-2.44391409, -3.40071697},
//      {-1.34133429,  0.18988270},
//      {-0.72723685,  2.57530980},
//      { 1.32475328,  3.89715240},
//      { 2.80434824,  3.23003402},
//      {-1.13351168,  0.62035750},
//      {-1.89567738,  0.37405691},
//      {-2.78032228,  0.69539255},
//      { 0.46263041,  3.13904618},
//      { 2.03036081,  4.98975253},
//      { 3.47932599,  7.47767211},
//      { 3.65911071,  7.27986112},
//      { 2.80631804,  7.34753688},
//      { 4.14584501,  8.61722424},
//      { 4.84461878,  5.35442857},
//      { 2.58185248,  1.98614034},
//      { 1.07744105,  0.23974312},
//      { 0.86648643, -2.11389155},
//      {-0.21132663, -2.66016188},
//      {-1.01162997, -1.40710159},
//      {-0.19728246, -1.18104782},
//      {-0.08026505, -1.95809150},
//      { 0.40098679, -1.30425030},
//      { 0.37458600, -0.45047353},
//      {-1.22499303,  0.06781530},
//      {-0.44528094,  0.68091861},
//      { 0.14750044,  2.00309904},
//      {-0.19018599,  3.01424746},
//      {-0.94157381,  3.93377398},
//      { 0.37093200,  8.01299089},
//      { 2.87473372, 10.32208970},
//      { 4.46626467,  9.20143323},
//      { 3.99070427,  7.71472151},
//      { 2.01383024,  5.58856868},
//      {-0.20770121,  4.88095156},
//      { 1.11814377,  5.74929748},
//      { 1.26080830,  5.17006810},
//      { 0.97129355,  6.04746852},
//      { 2.32421082,  6.33352591},
//      { 3.34670169,  4.77732191},
//      { 2.04275186,  3.25013277},
//      { 1.69112905,  3.32903825},
//      { 2.96856734,  2.65713531},
//      { 2.78884130,  0.99475339},
//      { 0.70796650, -1.31386473},
//      {-1.68837160, -0.39047899},
//      {-0.78381976,  0.52529795},
//      { 0.30155150,  0.45864152},
//      { 1.39062393,  0.91099202},
//      { 1.57032879,  0.52033224},
//      { 0.71873294, -1.34540456},
//      {-0.08451492, -2.66771179},
//      {-0.25904592, -6.03009777},
//      {-1.97302446, -6.37463119},
//      {-2.30676681, -5.62386338},
//      {-1.09087373, -6.07765782},
//      {-0.75975928, -4.70908817},
//      {-1.50542112, -1.13928850},
//      {-1.40117576,  1.96596256},
//      {-0.65233062,  2.68650764},
//      { 1.36185496,  2.97240640},
//      { 2.81503863,  1.55983671},
//      { 1.51706599, -0.98478712},
//      {-0.40306095, -2.06282504},
//      {-1.48787062, -1.84778663},
//      { 0.14580454, -2.02088053},
//      {-2.18155635, -2.18152051},
//      {-4.45305677,  0.10021396},
//      {-2.11506222,  5.09117423},
//      {-0.14560773,  7.28082157},
//      { 1.75277390,  9.34070122},
//      { 2.18443938, 11.52243600},
//      { 3.29005827, 12.41159116},
//      { 5.17717456, 10.78267370},
//      { 5.62433010,  8.69438371},
//      { 3.54458560,  7.25407699},
//      { 2.35591956,  7.00227126},
//      { 2.93497907,  6.02791002},
//      { 4.80080840,  3.94009985},
//      { 3.44469254,  1.64219866},
//      {-0.26822957, -0.00421737},
//      {-1.62324076,  0.27160381},
//      {-1.39268361,  0.71262399},
//      { 0.11080342,  3.65220609},
//      { 2.09347025,  5.77980333},
//      { 2.99392028,  5.42254057},
//      { 1.74440339,  5.05982814},
//      { 2.63642648,  7.43430643},
//      { 3.49151978,  7.94871742},
//      { 1.91464053,  7.94419568},
//      { 1.50544784,  8.78247786},
//      { 1.35179270, 10.42926215},
//      { 3.23317421, 11.68617769},
//      { 4.10497557, 10.37827843},
//      { 4.10834586,  8.75377476},
//      { 4.53683808,  8.53942825},
//      { 3.60077373,  6.94042807},
//      { 4.08295065,  4.31399431},
//      { 2.89419009,  1.44072132},
//      { 1.39565928, -0.89827551},
//      { 0.55197615, -1.07695549},
//      {-0.64374738, -2.41241436},
//      {-0.17868258, -1.86343908},
//      {-0.94382637, -1.21798914},
//      { 0.41745038, -1.15523075},
//      { 1.67111066, -2.80146958},
//      { 0.61357571, -3.64039574},
//      { 0.33363444, -4.51509156},
//      { 0.06105702, -3.67409972},
//      { 0.37736625, -4.07361907},
//      {-2.21637836, -6.29170828},
//      {-4.01018311, -6.24254618},
//      {-2.20576260, -4.68943349},
//      {-1.09687777, -3.86005608},
//      {-0.77072425, -3.46390869}};
//   CMatrix ts;
//   ts.elements = v_ts;
//   
//   CMatrix mask;
//   vector<double> para;
//   
//   int ar = 1;
//   int ma = 1;
//   bool ismean = false;
//   
//   Varma varma(ts, mask, para, ar, ma, ismean);
//   cout << "Hello World!";
//}


// MARK: RCpp Interfaces

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


RcppExport SEXP GetSVarmaResiduals(SEXP _timeSeries,
                                   SEXP _mask,
                                   SEXP _paramFixed,
                                   SEXP _orders,
                                   SEXP _arlags,
                                   SEXP _malags,
                                   SEXP _sresi,
                                   SEXP _swi,
                                   SEXP _isMeanIncluded
                                   ) {

   Rcpp::NumericMatrix RTimeSeries(_timeSeries);
   Rcpp::NumericMatrix RMask(_mask);
   Rcpp::NumericMatrix RSresi(_sresi);

   CMatrix TimeSeries = ToCMatrix(RTimeSeries);
   CMatrix Mask = ToCMatrix(RMask);
   CMatrix Sresi = ToCMatrix(RSresi);

   std::vector<int> Orders = Rcpp::as< std::vector<int> >(_orders);
   std::vector<int> ARLags = Rcpp::as< std::vector<int> >(_arlags);
   std::vector<int> MALags = Rcpp::as< std::vector<int> >(_malags);

   std::vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< std::vector<double> >(_paramFixed);
   }

   bool isMeanIncluded = Rcpp::as<bool>(_isMeanIncluded);
   bool swi = Rcpp::as<bool>(_swi);

   SVarma svarma(TimeSeries, Mask, ParamFixed, isMeanIncluded, Orders, ARLags, MALags, Sresi, swi);

   CMatrix Residuals = svarma.Residuals;
   return Rcpp::wrap(Residuals.elements);
}

RcppExport SEXP GetVMAObs(SEXP _timeSeries,
                          SEXP _mask,
                          SEXP _paramFixed,
                          SEXP _q,
                          SEXP _isMeanIncluded
                          ) {

   Rcpp::NumericMatrix RTimeSeries(_timeSeries);

   CMatrix TimeSeries = ToCMatrix(RTimeSeries);
   CMatrix Mask;
   if (!Rf_isNull(_mask)){
      Rcpp::NumericMatrix RMask(_mask);
      Mask = ToCMatrix(RMask);
   }

   std::vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< std::vector<double> >(_paramFixed);
   }

   int ma_q = Rcpp::as<int>(_q);

   bool isMeanIncluded = Rcpp::as<bool>(_isMeanIncluded);
   
   VMADemean VMADemean(TimeSeries, Mask, ParamFixed, ma_q, isMeanIncluded);
   return Rcpp::wrap(VMADemean.Obs.elements);
}

RcppExport SEXP GetVMATH(SEXP _timeSeries,
                         SEXP _mask,
                         SEXP _paramFixed,
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
   
   int ma_q = Rcpp::as<int>(_q);
   
   bool isMeanIncluded = Rcpp::as<bool>(_isMeanIncluded);
   
   VMA VMA(TimeSeries, Mask, ParamFixed, ma_q, isMeanIncluded);
   
   return Rcpp::wrap(VMA.TH.elements);
}
