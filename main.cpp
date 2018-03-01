// compile with
// c++ -o cm class_matrix.cpp -lcurl -std=c++11 -fopenmp
// execute with
// ./cm
#include"matrix.h"
#include<iostream>
#include<fstream>
#include<sstream>
#include<random>
#include<set>
#include<vector>
#include<algorithm>
#include<numeric>
#include<chrono>
#include<utility>
#include<cmath>
#include"hasting.h"

using std::exp;
using std::pair;
using std::partial_sum;
using std::vector;
using std::istringstream;
using std::string;
using std::set;
using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;



int main() {
	set<string> stocklist;
	stock_list(stocklist);
	print_all(stocklist);
	//clean_list(stocklist);
	cout<<"how many stocks that you must have"<<endl;
	int numofmust;
	cin>>numofmust;
	cout<<"what are they ?"<<endl;
	string ticker;
	set<string> must;
	for(int i=0;i<numofmust;){
		cin>>ticker;
		must.insert(ticker);
		auto flag=stocklist.find(ticker);
		if(flag==stocklist.end()){
			stockDataToFile(ticker);
			double* temp=read_data(ticker);
			if(temp[0]<pow(10,-4)&&temp[103]<pow(10,-4))
			{
				cout<<"can not get the proper data, please choose another one"<<endl;
				continue;
			}
			stocklist.insert(ticker);
		}
		++i;
	}
cout<<"How many stocks do you want to have in your portfolio ?"<<endl;
int num;
cin>>num;
cout<<"what is your expected return ?"<<endl;
int r;
cin>>r;
	//================daily data====================
	matrix data = matrix_maker(stocklist);
	//cout<<data<<endl;
//=======================calculate the covariance===
	matrix covariance = cov(data);
//===============calculate mean============
matrix m=mean(data);
int ssz=stocklist.size();
vector<string> ssname(ssz);
int ii=0;
for(auto i=stocklist.begin();i!=stocklist.end();++i){
	ssname[ii]=*i;
	++ii;
}
//============================================
pair<double, matrix> result=HM(covariance,m,stocklist, must, num,r);
cout<<"the lowest variance should be: "<<result.first<<"%"<<endl;
cout<<"the weight of this portfolio is: "<<endl;
matrix p=result.second;
int col=p.getcol();
for(int i=0;i<col;++i){
	cout<<ssname[p[0][i]]<<": "<<(p[1][i]*100.0)<<"%"<<endl;
}
	return 0;
}
