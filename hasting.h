#ifndef HASTING
#define HASTING
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

template<typename T>
void print_all(const T& v1){
	for(auto &c:v1){
		cout<<c<<" ";
	}
	cout<<endl;
}

void stock_list(set<string>& list) {
	ifstream in;
	string line;
	in.open("./quotes/stlist.txt");
	if (in.is_open()) {
		while (getline(in, line)) {
			list.insert(line);
		}
	}
	else {
		cout << "this file does not exist" << endl;
	}
	in.close();
}

double* read_data(const string& ticker) {
	ifstream in;
	string data;
	string line;
	double* result = new double[104];
	//vector<double> result;
  string address="./quotes/";
	address+=ticker;
	address+=".txt";
	in.open(address);
	int pos = 0;
	int count = 0;
	if (in.is_open()) {
		while (getline(in, line)) {
			++count;
			if (count == 1) {
				continue;
			}
			else {
				istringstream record(line);
				int tempc = 0;
				while (getline(record, data, ',')) {
					++tempc;
					if (tempc == 2) {
						//result.push_back(atof(data.c_str()));
						result[pos] = atof(data.c_str());
						++pos;
					}
				}
			}
		}

	}
	else {
		cout<<ticker<<endl;
		cout << "fail to open that file";
	}
	return result;
}

void clean_list(const set<string>& stocklist){
	ofstream out;
	out.open("./quotes/stlist.txt");
	for(auto& c:stocklist){
		double* result=read_data(c);
		if (result[0]>pow(10,-4)&&result[103]>pow(10,-4)){
			out<<c<<"\n";
		}
	}
	out.close();
}

matrix matrix_maker(set<string> &list) {
	int row = list.size();
	matrix data(row, 104);
	int i = 0;
	for (auto c = list.begin(); c != list.end(); ++c) {

		double* temp = read_data(*c);

		double* tt = data[i];
		data[i] = temp;
		temp = nullptr;
		delete[] tt;
		++i;
	}
	data=data*100.0;
	return std::move(data);
}

vector<int> get_index(const set<string> &stocklist, const set<string>& must)
{
	vector<int> mmust;
	auto startpoint = stocklist.begin();
	int index = -1;
	int pos = 0;
	for (auto j = must.begin(); j != must.end(); ++j)
	{
		for (auto i = startpoint; i != stocklist.end(); ++i) {
			++index;
			if (*j == *i) {
				mmust.push_back(index);
				i = ((++i) != stocklist.end()) ? i : --i;
				startpoint = i;
				break;
			}
		}
	}
	return std::move(mmust);
}

vector<int> get_other(const set<string>& stocklist, const vector<int>& mmust) {
	int numofuni = stocklist.size();
	int numofmust = mmust.size();
	vector<int>other;
	int numofignore = 0;
	for (int i = 0; i<numofuni; ++i) {
		int temp = 1;
		if (numofignore != numofmust) {
			for (auto j = mmust.begin(); j != mmust.end(); ++j) {
				if (i == *j)
				{
					temp = 0;
					++numofignore;
					break;
				}
			}
		}
		if (temp == 1) {
			other.push_back(i);
		}
	}
	return std::move(other);
}

pair<int, int> random_choose_stock(const vector<int>& other)
{
	int numofother = other.size();
	vector<double> tempprb(numofother, 1.0 / numofother);
	vector<double> prb(numofother);
	partial_sum(tempprb.begin(), tempprb.end(), prb.begin());
	//=============================generate seed==============================
	unsigned int seed2 = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine e2(seed2);
	std::uniform_real_distribution<double>u(0.0, 1.0);
	double p2 = u(e2);
	//============================================================================
	int stock;
	int pos = 0;
	for (int i = 1; i<numofother; ++i) {
		if (p2 < prb[0]) {
			stock = other[0];
			pos = 0;
			break;
		}
		else {
			if (p2 > prb[i - 1] && p2 <= prb[i])
			{
				stock = other[i];
				pos = i;
				break;
			}
		}
	}
	pair<int, int> result{ pos, stock };
	return std::move(result);
}

vector<int> choose_nbor(vector<int> mmust,  vector<int> ourchoice, const vector<int>& other) {
	int stock;
	do {
		stock = random_choose_stock(other).second;
	} while (find(ourchoice.begin(), ourchoice.end(), stock) != ourchoice.end());

	pair<int,int> takeout = random_choose_stock(ourchoice);
	ourchoice[takeout.first] = stock;

	return std::move(ourchoice);
}


pair<double, matrix> min_cov(matrix& cov, matrix& mean, vector<int> mmust, vector<int> ourchoice, const double& miu)
{
	/*cout<<"I AM NOW IN THE FUCNTION MIN_COV"<<endl;
	cout<<"the stock that user must have are: "<<endl;
	print_all(mmust);
	cout<<"ourchoice"<<endl;
	print_all(ourchoice);*/
	mmust.insert(mmust.begin(), ourchoice.begin(), ourchoice.end());
	sort(mmust.begin(), mmust.end());
	/*cout<<"we already have: "<<endl;
	print_all(mmust);*/
	//==========================================
	int sz = mmust.size();
	matrix sub_cov(sz, sz);
	matrix sub_mean(1, sz);
	for (int i = 0; i<sz; ++i) {
		int pos1 = mmust[i];
		sub_mean[0][i] = mean[pos1][0];
		for (int j = 0; j <= i; ++j) {
			int pos2 = mmust[j];
			sub_cov[i][j] = cov[pos1][pos2];
			sub_cov[j][i] = cov[pos2][pos1];
		}
	}

	matrix A(1, 2);
	A[0][0] = miu;
	A[0][1] = 1;
	/*cout<<"A is: "<<endl;
	cout<<A<<endl;*/

	matrix invC = inv(sub_cov);
	//cout<<"invC is: "<<endl;
	//cout<<invC<<endl;

	matrix sub_meant = sub_mean;
	sub_meant.T();
	/*cout<<sub_mean.getrow()<<" "<<sub_mean.getcol()<<endl;
	cout<<invC.getrow()<<" "<<invC.getcol()<<endl;*/
	matrix tempB = sub_mean*invC*sub_meant;
	double B = tempB[0][0];
	matrix u(1, sz, 1);
	matrix ut = u;
	ut.T();
	matrix tempC = sub_mean*invC*ut;
	double C = tempC[0][0];
	matrix tempD = u*invC*sub_meant;
	double D = tempD[0][0];
	matrix tempE = u*invC*ut;
	double E = tempE[0][0];
	matrix M(2, 2);
	M[0][0] = B; M[0][1] = C; M[1][0] = D; M[1][1] = E;
	matrix invM = inv(M);
	matrix ttemp = A*invM;
	matrix tempF = sub_mean*invC;
	matrix tempG = u*invC;
	matrix wgt = tempF*ttemp[0][0] + tempG*ttemp[0][1];
	matrix wgtT = wgt;
	wgtT.T();
	matrix tempcov = wgt*sub_cov*wgtT;
	double mincov = tempcov[0][0];
	pair<double, matrix>result{ mincov, wgt };
	return std::move(result);
}

pair<double, matrix> HM(matrix& cov,matrix& mean, set<string> &stocklist, set<string>& must, const int& num, const double& miu){
	int numofuni=stocklist.size();
	int numofmust=must.size();
	//===================================================
	vector<int> mustindex=get_index(stocklist,must);
	//cout<<"the index of the stock that user must have is: "<<endl;
	//print_all(mustindex);
	//====================================================
	int numofourchoice=num-numofmust;
	vector<int> other = get_other(stocklist, mustindex);
	//cout<<"what else can we choose: "<<endl;
	//print_all(other);
	int numofother=other.size();
	//====================================================
	vector<int> ourchoice(numofourchoice);
	//=====================================================
	vector<int> can_choose=other;
	//================================================
	for(int i=0;i<numofourchoice;++i){
		//cout<<"what else can we choose: "<<endl;
		//print_all(can_choose);
		int stock=random_choose_stock(can_choose).second;
		remove(can_choose.begin(), can_choose.end(), stock);
		can_choose.pop_back();
		ourchoice[i]=stock;
	}
	//cout<<"what do we choose at the first time: "<<endl;
	//print_all(ourchoice);
	//===============================================================
	double T=10.0;
	double N=numofourchoice*(numofuni-numofmust-numofourchoice);
	pair<double, matrix> lowerst_state=min_cov(cov,mean,mustindex,ourchoice, miu);
	//cout<<"we finished to calculate the lowerst_state"<<endl;
	double lowerst_energy=lowerst_state.first;
	//cout<<"first time lowerst_energy is: "<<lowerst_energy<<endl;
	double energy;
	double deltaE;
	double last_energy=lowerst_energy;
	pair<double, matrix> new_state;
	pair<double, matrix>last_state;
	vector<int>lowerst_choice;
	//(matrix& cov, matrix& mean, vector<int> mmust, vector<int>& ourchoice, const double& miu)
	int trytime=10000;
	for(int i=0;i<trytime;++i){
		//choose_nbor(vector<int> mmust,  vector<int> ourchoice, const vector<int>& other) {
		//cout<<"ourchoice： "<<endl;
		//print_all(ourchoice);
		vector<int>newourchoice=choose_nbor(mustindex, ourchoice, other);
		//cout<<"my new choose is: "<<endl;
		//print_all(newourchoice);
		new_state=min_cov(cov, mean, mustindex, newourchoice, miu);
		energy=new_state.first;
		//cout<<i+1<<", energy is: "<<energy<<endl;
		deltaE=energy-last_energy;
		//cout<<"deltaE is: "<<deltaE<<endl;
		if(deltaE<=0){
			//cout<<"deltaE "<<endl;
			//cout<<"update information"<<endl;
			ourchoice=std::move(newourchoice);
			//cout<<"ourchoice : "<<endl;
			//print_all(ourchoice);

			last_state=std::move(new_state);
			//cout<<"last_state becomes: "<<endl;
			//cout<<last_state.first<<endl;
			//cout<<last_state.second<<endl;
			last_energy=last_state.first;
			if(last_energy<lowerst_energy)
			{
				//cout<<"energy比lowerst energy"<<endl;
				lowerst_choice=ourchoice;
				lowerst_state=last_state;
				lowerst_energy=last_energy;
				//cout<<"so, lowerst_energy becomes: "<<endl;
				//cout<<lowerst_energy<<endl;
			}
		}
		else{
			//cout<<"deltaE >0 at this time"<<endl;
			double temp=(-deltaE)/T;
			double threshhold=exp(temp);
			//cout<<"threshhold is: "<<endl;
			//cout<<threshhold<<endl;
			unsigned int seed2 = std::chrono::system_clock::now().time_since_epoch().count();
			std::default_random_engine e2(seed2);
			std::uniform_real_distribution<double>u(0.0, 1.0);
			double p2 = u(e2);
			//cout<<"the random number we generate is: "<<p2<<endl;
			if(p2<threshhold){
				//cout<<"p2<0, state"<<endl;
				ourchoice=std::move(newourchoice);
				//cout<<"ourchoice : "<<endl;
				//print_all(ourchoice);
				last_state=std::move(new_state);
				last_energy=last_state.first;
				//cout<<"last_energy : "<<endl;
				//cout<<last_energy<<endl;
			}
		}
	}
	mustindex.insert(mustindex.begin(), lowerst_choice.begin(), lowerst_choice.end());
	sort(mustindex.begin(), mustindex.end());
	matrix result(2,num);
	double* temp=result[1];
	result[1]=lowerst_state.second[0];
	lowerst_state.second[0]=nullptr;
	delete[] temp;

	for(int i=0;i<num;++i){
		result[0][i]=mustindex[i];
	}
	pair<double, matrix> finalresult{lowerst_state.first,result};
	return std::move(finalresult);
}



#endif
