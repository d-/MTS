#include <vector>
#include "Rcpp.h"

using namespace std;
using namespace Rcpp;

// MARK: CMatrix
class CMatrix {
public:
   vector< vector<double> > elements;
   
   CMatrix(double ele, int n_row, int n_col);
   CMatrix(vector< vector<double> > eles);
   CMatrix();
   
   int nrow();
   int ncol();
   
   double * operator()(int p, int q);
   vector<double> operator() (int p, bool is_row);
   
   void append(vector<double>);
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

CMatrix::CMatrix(vector< vector<double> > eles) {
   elements = eles;
}

CMatrix::CMatrix(double ele, int n_row, int n_col) {
   for (int i = 0; i < n_row; i++){
      vector<double> temp_col(n_col);
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
      cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
}

void CMatrix::append(vector<double> row){
   if (ncol() == 0 || ncol() == row.size()){
      elements.push_back(row);
   } else {
      cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
   }
}

void CMatrix::transpose() {
   vector< vector<double> > transposed_elements;
   if (elements.size() > 0){
      for (int j = 0; j < elements[0].size(); j++) {
         vector<double> temp;
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

vector<double> CMatrix::operator()(int p, bool is_row = true) {
   vector<double> ret;
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
         vector<double> temp;
         for (int i = 0; i < A.elements.size(); i++) {
            temp.push_back(A.elements[i][j]);
         }
         B.elements.push_back(temp);
      }
   }
   return (B);
}

CMatrix Cdiagonal (vector<double> &A) {
   CMatrix B (0, (int) A.size(), (int) A.size());
   for (int i = 0; i < A.size(); i++){
      B.elements[i][i] = A[i];
   }
   return (B);
}

CMatrix Cidentity (int s){
   vector<double> ones(s);
   fill(ones.begin(), ones.end(), 1);
   CMatrix A = Cdiagonal(ones);
   return (A);
}

CMatrix ToCMatrix (Rcpp::NumericMatrix RMat){
   CMatrix C;
   int nrow = RMat.nrow();
   int ncol = RMat.ncol();
   for (int i = 0; i < nrow; i++){
      vector<double> row;
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
      cerr<<"rbind Error:  Can't rbind matrices with different column sizes.\n";
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
      cerr<<"rbind Error:  Can't cbind matrices with different row sizes.\n";
   }
   
   return C;
}



CMatrix as_matrix (vector<double> x, bool is_as_row = true){
   CMatrix C;
   if(is_as_row){
      C.elements.push_back(x);
   } else {
      for (int i = 0; i < x.size(); i++){
         vector<double> t;
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
      cerr<<"CMatrix Prod Error:  Incompatible matrices.\n";
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
   
   vector<double> Ph0; // Constants
   CMatrix PH; // AR Coeff
   CMatrix TH; // MA Coeff
   
   int p;  // AR order p
   int q;  // MA order q
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() ==0
   
   Varma(CMatrix &, CMatrix &, vector<double> &, int, int, bool);
   
private:
   int checkMaskFormat(CMatrix &);
   void fillParamFixed(CMatrix &, vector<double>, bool);
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
            cerr<<"Invalid Mask CMatrix:  CMatrix contains elements others than 0 or 1 ("
            <<i<<", "<<j<<") = "<<*Mask(i,j)<<"\n";
            break;
         }
      }
   }
   return sum;
}

void Varma::fillParamFixed(CMatrix & Mask, vector<double> ParamFixed, bool isMeanIncluded) {
   CMatrix Beta;
   int kp = k * p;
   int kq = k * q;
   
   int i_start = 0;
   
   vector<double> QParamFixed(ParamFixed.size());
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
      cerr<<"Init with parameters error:  Too many parameters\n";
   }
   
}

void Varma::compResiduals(){
   // Const Row
   vector<double> Res_Const;
   for (int i = 1; i <= k; i++){
      Res_Const.push_back(*Obs(1,i)-Ph0[i-1]);
   }
   
   Residuals.append(Res_Const);
   
   // ------------------------------------------------------------------------
   // Step 1:
   // Calculate max(p,q) residuals, necessary for subsequent recursions
   if (max(p,q)>1){
      for (int t = 2; t<= max(p,q); t++){
         vector<double> Res_Row;
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
               
               // Downcasting matrix to row vector
               vector<double> Estimate_AR = prod(as_matrix(Obs(t-j)),PH_Slice)(1);
               
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
               
               // Downcasting matrix to row vector
               vector<double> Estimate_MA = prod(as_matrix(Residuals(t-j)),TH_Slice)(1);
               
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
   
   int i_start = max(p,q)+1;
   for (int t = i_start; t<= nT; t++) {
      vector<double> Obs_Past;
      if (p > 0) {
         for (int j=1; j<=p; j++){
            vector<double> Zt_Slice = Obs(t-j);
            Obs_Past.insert(Obs_Past.end(), Zt_Slice.begin(),Zt_Slice.end());
         }
      }
      
      if (q > 0) {
         for (int j=1; j<=q; j++){
            vector<double> At_Slice = Residuals(t-j);
            Obs_Past.insert(Obs_Past.end(), At_Slice.begin(),At_Slice.end());
         }
      }
      
      Obs_Past.insert(Obs_Past.begin(), Obs_Const);
      
      vector<double> Estimate_ARMA = prod(as_matrix(Obs_Past), Beta)(1);
      
      vector<double> Res_Row;
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
Varma::Varma (CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed, int ar_p, int ma_q, bool isMeanIncluded) {
   Obs = TimeSeries;
   
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   
   p = ar_p;
   q = ma_q;
   
   // Check compactible:  Mask and Parameters
   int mask_size = checkMaskFormat(Mask);
   
   if (mask_size != ParamFixed.size()) {
      cerr<<"Input error:  Fixed mask and parameters provided are of different size.\n";
      cerr<<"   Mask contains "<<mask_size<<" 1s but "<<ParamFixed.size()<<" is provided.";
   } else {
      // Check compactible:  Mask and p, q, k
      int row_const = hasMean? 1 : 0;
      
      if (row_const+k*(p+q) != Mask.nrow()) {
         cerr<<"Input error:  Fixed mask has a wrong row number.\n";
         cerr<<"Need "<<row_const+k*(p+q)<<" rows but the mask matrix has "<<Mask.nrow()<<" rows.\n";
      } else {
         fillParamFixed(Mask, ParamFixed, hasMean);
         
         compResiduals();
      }
   }
}

// MARK: SVarma
class SVarma {
public:
   CMatrix Obs;
   CMatrix Residuals;
   int k;  // # Param
   int nT; // # Obs
   
   vector<double> Ph0; // Constants
   CMatrix PH; // AR Coeff
   CMatrix TH; // MA Coeff
   CMatrix sPH; // sAR Coeff
   CMatrix sTH; // sMA Coeff
   
   bool matrix_prod_method_switch; // swi
   CMatrix Phi; //
   CMatrix Theta; //

   CMatrix Beta; //!!! Local property in Varma
   vector<int> ARlags;
   vector<int> MAlags;
   
   int nar;  // # Param
   int nma; // # Obs
   
   int p;  // AR order p
   int q;  // MA order q
   int P;  // SAR order P
   int Q;  // SMA order Q
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() ==0
   
   SVarma (CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed, bool isMeanIncluded, vector<int> Orders, vector<int> ARlags, vector<int> MAlags, CMatrix & Sresi, bool swi);
private:
   int checkMaskFormat(CMatrix & Mask);
   void fillParamFixed(CMatrix & Mask, vector<double> ParamFixed, bool isMeanIncluded);
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
            cerr<<"Invalid Mask CMatrix:  CMatrix contains elements others than 0 or 1 ("
            <<i<<", "<<j<<") = "<<*Mask(i,j)<<"\n";
            break;
         }
      }
   }
   return sum;
}

void SVarma::fillParamFixed(CMatrix & Mask, vector<double> ParamFixed, bool isMeanIncluded) {
   CMatrix Beta;
   int kp = k * p;
   int kq = k * q;
   int kP = k * P;
   int kQ = k * Q;
   
   int i_start = 0;
   
   vector<double> QParamFixed(ParamFixed.size());
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
         cerr<<"Init with parameters error:  Too many parameters\n";
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
         cerr<<"Init with parameters error:  Too many parameters\n";
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

SVarma::SVarma (CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed, bool isMeanIncluded, vector<int> Orders, vector<int> ARlags, vector<int> MAlags, CMatrix & Sresi, bool swi)
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
   i_start = max(ARlags.back(), MAlags.back())+1;
   
   
   // Check compactible:  Mask and Parameters
   int mask_size = checkMaskFormat(Mask);
   
   if (mask_size != ParamFixed.size()) {
      cerr<<"Input error:  Fixed mask and parameters provided are of different size.\n";
      cerr<<"   Mask contains "<<mask_size<<" 1s but "<<ParamFixed.size()<<" is provided.";
   } else {
      // Check compactible:  Mask and p, q, k
      int row_const = hasMean? 1 : 0;
      
      if (row_const+k*(nar+nma) != Mask.nrow()) {
         cerr<<"Input error:  Fixed mask has a wrong row number.\n";
         cerr<<"Need "<<row_const+k*(p+q)<<" rows but the mask matrix has "<<Mask.nrow()<<" rows.\n";
      } else {
         fillParamFixed(Mask, ParamFixed, hasMean);
         compResiduals(Sresi);
      }
   }
   
   
}


// MARK: VMA
class VMA {
public:
   CMatrix Obs;
   CMatrix Residuals;
   int k;  // # Param
   int nT; // # Obs
   
   int q;  // MA order q

   vector<double> Ph0; // Constants -- OC: mu

   CMatrix Theta;
   CMatrix TH;

   bool hasMean;  // Whether mean is included, i.e. Ph0.size() == 0
   
   VMA(CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed,int ma_q, bool isMeanIncluded);
private:
};

VMA::VMA(CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed,int ma_q, bool isMeanIncluded)
{
   Obs = TimeSeries;
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   q = ma_q;
   
   vector<double> QParamFixed(ParamFixed.size());
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
   
   vector<double> Ph0; // Constants -- OC: mu
   
   CMatrix Theta;
   CMatrix TH;
   
   bool hasMean;  // Whether mean is included, i.e. Ph0.size() == 0
   
   VMADemean(CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed,int ma_q, bool isMeanIncluded);
private:
};

VMADemean::VMADemean(CMatrix & TimeSeries, CMatrix & Mask, vector<double> & ParamFixed,int ma_q, bool isMeanIncluded)
{
   Obs = TimeSeries;
   k = (int) Obs.ncol();
   nT = (int) Obs.nrow();
   hasMean = isMeanIncluded;
   q = ma_q;
   
   vector<double> QParamFixed(ParamFixed.size());
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




// MARK: Interfaces

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
   
   vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< vector<double> >(_paramFixed);
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
   
   vector<int> Orders = Rcpp::as< vector<int> >(_orders);
   vector<int> ARLags = Rcpp::as< vector<int> >(_arlags);
   vector<int> MALags = Rcpp::as< vector<int> >(_malags);

   vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< vector<double> >(_paramFixed);
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
   Rcpp::NumericMatrix RMask(_mask);
   
   CMatrix TimeSeries = ToCMatrix(RTimeSeries);
   CMatrix Mask = ToCMatrix(RMask);
   
   vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< vector<double> >(_paramFixed);
   }
   
   int ma_q = Rcpp::as<int>(_q);

   bool isMeanIncluded = Rcpp::as<bool>(_isMeanIncluded);
   
   VMADemean VMADemean(TimeSeries, Mask, ParamFixed, ma_q, isMeanIncluded);
   return wrap(VMADemean.Obs.elements);
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
   
   vector<double> ParamFixed;
   if (!Rf_isNull(_paramFixed)){
      ParamFixed= Rcpp::as< vector<double> >(_paramFixed);
   }
   
   int ma_q = Rcpp::as<int>(_q);
   
   bool isMeanIncluded = Rcpp::as<bool>(_isMeanIncluded);
   
   VMA VMA(TimeSeries, Mask, ParamFixed, ma_q, isMeanIncluded);
   
   return wrap(VMA.TH.elements);
}

