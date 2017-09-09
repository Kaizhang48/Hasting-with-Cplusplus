#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include<vector>
#include "omp.h"
#include <string>
#include <fstream>
#include <curl/curl.h>
using std::vector;
using std::ostream;
using std::istream;
using std::cout;
using std::cin;
using std::endl;
using std::out_of_range;

typedef long int myint;
void printToFile(std::string filename, const std::string &text){
std::ofstream mfile;
mfile.open(filename);
for(myint i=0;i<text.length();++i){
mfile<<text[i];
}
mfile.close();
}
static size_t writerF(void *ptr, size_t size, size_t nmemb, void *userdata)
{
((std::string*)userdata)->append((char*)ptr, size * nmemb);
return size * nmemb;
}
void stockDataToFile(const std::string &tickerName,
const std::string &quandl_auth_token="d7BCxzFbxDeVAPSgYjwi",
const std::string &database="WIKI",
const std::string &folder="./quotes/"){
/*std::string mainLink="https://www.quandl.com/api/v1/datasets/";
mainLink+="/"+tickerName;
mainLink+=".csv";
mainLink+="?sort_order=asc&auth_token=";
*/
std::string mainLink="https://www.quandl.com/api/v3/datasets/";
mainLink+=database;
mainLink+="/"+tickerName;
mainLink+=".csv?column_index=4&start_date=2016-12-05&end_date=2017-05-05&collapse=daily&transform=rdiff&api_key=";
mainLink+=quandl_auth_token;
CURL *curl;
std::string quandlData;
std::string fName=folder;
fName+=tickerName;
fName+=".txt";
curl = curl_easy_init();
if(curl) {
const char* linkArrayChar=mainLink.c_str();
curl_easy_setopt(curl, CURLOPT_URL, linkArrayChar);
curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writerF);
curl_easy_setopt(curl, CURLOPT_WRITEDATA, &quandlData);
curl_easy_perform(curl);
curl_easy_cleanup(curl);
printToFile(fName,quandlData);
}
}

class matrix {
	friend matrix mean(const matrix& );
	friend matrix inv(matrix a);
	friend matrix operator* (const matrix& A,const matrix& B);
	friend matrix operator* (const matrix& A, const double &B);
	friend matrix operator+ (const matrix& A,  const matrix &B);
	friend matrix operator- (const matrix& A,  const matrix& B);
	friend ostream& operator <<(std::ostream &os, const matrix &m);
public:
	double*& operator[](int t);
	matrix() = default;
	matrix(double* const a,const int&n);
	matrix(const int& a, const int& b,const double& c=0);
	matrix(const matrix &copyfrom);
	matrix(matrix &&movefrom);
	~matrix();
	int getrow()const;
	int getcol()const;
	double getdata(const int&i,const int&j)const;
	matrix& operator=(const matrix& assignfrom);

	matrix& operator=(matrix&& moveassignfrom);
	vector<int> size() const ;
	matrix& T();

private:
	int row;
	int col;
	double** data;
	void mfree();
};

void matrix::mfree() {
  if (data != nullptr) {
#pragma omp parallel
{
#pragma omp for
    for (int i = 0; i < row; ++i)
		{
      delete[] data[i];
		}
}
    delete[] data;
  }
  row = 0;
  col = 0;
}

double*& matrix::operator[](int t) {
  double* &r = data[t];
  return r;
}
matrix::matrix(double* const a,const int&n):row(1),col(n) {
  data = new double*[1];
  data[0] = new double[n];
  for(int i=0;i<col;++i){
    data[0][i]=a[i];
  }
}
matrix::matrix(const int&a, const int& b,const double& c) :row(a), col(b) {
  data = new double*[row];
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < row; ++i)
      data[i] = new double[col];
  }
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < row; ++i)
    {
      for (int j = 0; j < col; ++j)
      {
        data[i][j] = c;
      }
    }
  }
}
matrix::matrix(const matrix &copyfrom)
{
  //cout << "this is copy constructor" <<endl;
  row = copyfrom.row;
  col = copyfrom.col;
  data = new double*[row];
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < row; ++i) {
      //printf("I am Thread %d\n", omp_get_thread_num());
      data[i] = new double[col];
    }
  }
#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < row; ++i)
    {
      //printf("I am Thread %d\n", omp_get_thread_num());
      for (int j = 0; j < col; ++j)
      {
        data[i][j] = copyfrom.data[i][j];
      }
    }
  }
}
matrix::matrix(matrix &&movefrom) {
  if (movefrom.data != data) {
    data = movefrom.data;
    row = movefrom.row;
    col = movefrom.col;
    movefrom.data = nullptr;
    movefrom.row = 0;
    movefrom.col = 0;
  }
}
matrix& matrix::operator=(const matrix& assignfrom) {
  if (assignfrom.data != data) {
    if (data != nullptr) {
      mfree();
    }
    row = assignfrom.row;
    col = assignfrom.col;
    data = new double*[row];
#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < row; ++i)
        data[i] = new double[col];
    }
#pragma omp parallel
    {
#pragma omp for
      for (int i = 0; i < row; ++i)
      {
        for (int j = 0; j < col; ++j)
        {
          data[i][j] = assignfrom.data[i][j];
        }
      }
    }
  }
	return *this;
}

matrix& matrix::operator=(matrix&& moveassignfrom) {
  if (moveassignfrom.data != data) {
    if (data != nullptr) {
      mfree();
    }
    data = moveassignfrom.data;
    row = moveassignfrom.row;
    col = moveassignfrom.col;
    moveassignfrom.data = nullptr;
    moveassignfrom.row = 0;
    moveassignfrom.col = 0;
  }
	return *this;
}

matrix::~matrix()
{
   mfree();
 }

 int matrix::getrow()const{
   return row;
 }

 int matrix::getcol()const{
   return col;
 }

 double matrix::getdata(const int&i,const int&j)const{
   return data[i][j];
 }

 vector<int> matrix::size() const {
   return vector<int>{row, col};
 }

 matrix& matrix::T() {
   double** tempdata = new double*[row];
#pragma omp parallel
   {
#pragma omp for
     for (int i = 0; i < row; ++i) {
       //printf("I am Thread %d\n", omp_get_thread_num());
       tempdata[i] = new double[col];
     }
   }
#pragma omp parallel
   {
#pragma omp for
     for (int i = 0; i < row; ++i)
     {
       //printf("I am Thread %d\n", omp_get_thread_num());
       for (int j = 0; j < col; ++j)
       {
         tempdata[i][j] = data[i][j];
       }
     }
   }
   int tempr = col;
   int tempc = row;
   mfree();
   row = tempr;
   col = tempc;
   data = new double*[row];
#pragma omp parallel
   {
#pragma omp for
     for (int i = 0; i < row; ++i) {
       //printf("I am Thread %d\n", omp_get_thread_num());
       data[i] = new double[col];
     }
   }

#pragma omp parallel
   {
#pragma omp for
     for (int i = 0; i < row; ++i) {
       for (int j = 0; j < col; ++j) {
         data[i][j] = tempdata[j][i];
       }
     }
   }
#pragma omp parallel
   {
#pragma omp for
     for (int i = 0; i < col; ++i) {
       delete[] tempdata[i];
     }
   }
   delete[] tempdata;
   return *this;
 }
 /*matrix inv(matrix b)
 {
 	double**&a = b.data;
 	int& n = b.row;
 	int *is = new int[n];
 	int *js = new int[n];
 	int i, j, k;
 	double d, p;
 	for (k = 0; k < n; k++)
 	{
 		d = 0.0;
 		for (i = k; i <= n - 1; i++)
 			for (j = k; j <= n - 1; j++)
 			{
 				p = fabs(a[i][j]);
 				if (p>d) { d = p; is[k] = i; js[k] = j; }
 			}
 		if (0.0 == d)
 		{
 			free(is); free(js);
 			throw out_of_range("can not be inversed !");
 		}
 		if (is[k] != k)
 			for (j = 0; j <= n - 1; j++)
 			{
 				p = a[k][j];
 				a[k][j] = a[is[k]][j];
 				a[is[k]][j] = p;
 			}
 		if (js[k] != k)
 			for (i = 0; i <= n - 1; i++)
 			{
 				p = a[i][k];
 				a[i][k] = a[i][js[k]];
 				a[i][js[k]] = p;
 			}
 		a[k][k] = 1.0 / a[k][k];
 		for (j = 0; j <= n - 1; j++)
 			if (j != k)
 			{
 				a[k][j] *= a[k][k];
 			}
 		for (i = 0; i <= n - 1; i++)
 			if (i != k)
 				for (j = 0; j <= n - 1; j++)
 					if (j != k)
 					{
 						a[i][j] -= a[i][k] * a[k][j];
 					}
 		for (i = 0; i <= n - 1; i++)
 			if (i != k)
 			{
 				a[i][k] = -a[i][k] * a[k][k];
 			}
 	}
 #pragma omp parallel
 {
 #pragma omp for
   	for (k = n - 1; k >= 0; k--)
   	{
   		if (js[k] != k)
   			for (j = 0; j <= n - 1; j++)
   			{
 #pragma omp critical
   				p = a[k][j];
   				a[k][j] = a[js[k]][j];
   				a[js[k]][j] = p;
   			}
   		if (is[k] != k)
   			for (i = 0; i <= n - 1; i++)
   			{
 #pragma omp critical
   				p = a[i][k];
   				a[i][k] = a[i][is[k]];
   				a[i][is[k]] = p;
   			}
   	}
 }
 	free(is); free(js);
 	return std::move(b);
 }
*/
matrix inv(matrix b)
{
 double**&a = b.data;
 int& n = b.row;
 int *is = new int[n];
 int *js = new int[n];
 int i, j, k;
 double d, p;
 for (k = 0; k < n; k++)
 {
	 d = 0.0;
	 for (i = k; i <= n - 1; i++)
		 for (j = k; j <= n - 1; j++)
		 {
			 p = fabs(a[i][j]);
			 if (p>d) { d = p; is[k] = i; js[k] = j; }
		 }
	 if (0.0 == d)
	 {
		 free(is); free(js);
		 throw out_of_range("can not be inversed !");
	 }
	 if (is[k] != k)
		 for (j = 0; j <= n - 1; j++)
		 {
			 p = a[k][j];
			 a[k][j] = a[is[k]][j];
			 a[is[k]][j] = p;
		 }
	 if (js[k] != k)
		 for (i = 0; i <= n - 1; i++)
		 {
			 p = a[i][k];
			 a[i][k] = a[i][js[k]];
			 a[i][js[k]] = p;
		 }
	 a[k][k] = 1.0 / a[k][k];
	 for (j = 0; j <= n - 1; j++)
		 if (j != k)
		 {
			 a[k][j] *= a[k][k];
		 }
	 for (i = 0; i <= n - 1; i++)
		 if (i != k)
			 for (j = 0; j <= n - 1; j++)
				 if (j != k)
				 {
					 a[i][j] -= a[i][k] * a[k][j];
				 }
	 for (i = 0; i <= n - 1; i++)
		 if (i != k)
		 {
			 a[i][k] = -a[i][k] * a[k][k];
		 }
 }


	 for (k = n - 1; k >= 0; k--)
	 {
		 if (js[k] != k)
			 for (j = 0; j <= n - 1; j++)
			 {
				 p = a[k][j];
				 a[k][j] = a[js[k]][j];
				 a[js[k]][j] = p;
			 }
		 if (is[k] != k)
			 for (i = 0; i <= n - 1; i++)
			 {
				 p = a[i][k];
				 a[i][k] = a[i][is[k]];
				 a[i][is[k]] = p;
			 }
	 }

 free(is); free(js);
 return std::move(b);
}

 matrix operator*( const matrix& A,const matrix& B) {
 	auto szA = A.size();
 	auto szB = B.size();
 	if (szA[1] != szB[0]) {
 		throw out_of_range("size of two matrixs does not match!");
 	}
 	matrix result(szA[0], szB[1]);
 	if(szA[0]==1){
	#pragma omp parallel
	{
	#pragma omp for
 		for(int i=0;i<szB[1];++i){
			for(int j=0;j<szA[1];++j){
	 	  	result.data[0][i]+=(A.getdata(0,j)*B.getdata(j,i));
			}
 		}
	}
 	}
 	else{
 #pragma omp parallel
 {
 #pragma omp for

 		for (int i = 0; i < szA[0]; ++i)
 		{
 			//printf("I am Thread %d\n", omp_get_thread_num());
 			for (int j = 0; j < szB[1]; ++j)
 			{
 				for (int k = 0; k < szA[1]; k++)
 				{
 #pragma omp critical
 					result.data[i][j] += (A.getdata(i,k)* B.getdata(k,j));
 				}
 			}
 		}
 }
 }
 	return std::move(result);
}

 matrix operator*( matrix& A, const double &B) {
 	auto sz = A.size();
 	matrix result = A;
 #pragma omp parallel
 	{
 #pragma omp for
 		for (int i = 0; i < sz[0]; ++i) {
 			for (int j = 0; j < sz[1]; ++j) {
 				result[i][j] *= B;
 			}
 		}
 	}
 	return std::move(result);
 }

 matrix operator*( const double &B, matrix& A){
 	return std::move(A*B);
 }

 matrix operator+ ( const matrix& A,  const matrix &B) {
 	auto szA = A.size();
 	auto szB = B.size();
 	if (szA[0] != szB[0] || szA[1] != szB[1]) {
 		throw out_of_range("size of two matrixs does not match!");
 	}
 	matrix c = A;
 #pragma omp parallel
 	{
 #pragma omp for
 			for (int i = 0; i < szA[0]; ++i) {
 				//printf("I am Thread %d\n", omp_get_thread_num());
 				for (int j = 0; j < szA[1]; ++j) {
 					c[i][j] += B.getdata(i,j);
 				}
 			}
 	}
 	return std::move(c);
 }
 matrix operator- ( matrix& A,  matrix& B) {
 	auto szA = A.size();
 	auto szB = B.size();
 	if (szA[0] != szB[0] || szA[1] != szB[1]) {
 		throw out_of_range("size of two matrixs does not match!");
 	}
 	matrix c = A;
 #pragma omp parallel
 	{
 #pragma omp for
 		for (int i = 0; i < szA[0]; ++i) {
 			//printf("I am Thread %d\n", omp_get_thread_num());
 			for (int j = 0; j < szA[1]; ++j) {
 				c[i][j] -= B.getdata(i,j);
 			}
 		}
 	}
 	return std::move(c);
 }

 ostream& operator <<(std::ostream &os, const matrix &m)
 {
 	for (int i = 0; i < m.row; i++)
 	{
 		os << " | ";
 		for (int j = 0; j < m.col; j++)
 		{
 			char buf[32];
 			double data = m.data[i][j];
 			if (m.data[i][j] > -0.00001 && m.data[i][j] < 0.00001)
 				data = 0;
 			sprintf(buf, "%10.10lf ", data);
 			os << buf;

 		}
 		os << "|\n";
 	}
 	os << "\n\n";
 	return os;
 }

 matrix mean(const matrix& a) {
	const vector<int> sz=a.size();
 	matrix result(sz[0], 1);
 #pragma omp parallel
 	{
 #pragma omp for
 		for (int i = 0; i < sz[0]; ++i) {
 			for (int j = 0; j < sz[1]; ++j) {
 				result[i][0] += a.data[i][j];
 			}
 			result[i][0]/=sz[1];
 		}
 	}
 	return std::move(result);
 }


 void removemean(matrix& a){
 	matrix b=mean(a);
 #pragma omp parallel
 {
 #pragma omp for
 	for (int i=0;i<a.getrow();++i){
 		for(int j=0;j<a.getcol();++j){
 			a[i][j]-=b[i][0];
 		}
 	}
 }
 }

matrix cov(matrix a){
		removemean(a);
		matrix aT=a;
		aT.T();
		matrix result=a*aT;
		result=result*(1.0/(a.getcol()-1));
		return std::move(result);
}
 /*matrix cov(matrix a){
 	removemean(a);
 	matrix result(a.getrow(), a.getrow());
 #pragma omp parallel
 {
 #pragma omp for
 	for (int i = 0; i < a.getrow(); ++i) {
 		for (int j = 0; j <= i; ++j) {
 			matrix tempi(a[i], a.getcol());
 			matrix tempj(a[j], a.getcol());
 			result[i][j] = result[j][i] = (tempi*tempj.T())[0][0]/(a.getcol()-1);
 		}
 	}
 }
 	return std::move(result);
}*/

#endif
