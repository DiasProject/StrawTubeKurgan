#include <iostream>
#include <vector>
#include "ctime"
using namespace std;
int main()
{
	srand(time(NULL));
    // Вектор из 10 элементов типа int
    vector<int> v1(10, 12);
	vector<double>	v4(12);
	for(int i=0; i<v1.size(); i++){
    cout<<"vector="<<v1[i]<<endl;
	}for(int i=0; i<v4.size(); i++){
	v4[i]=rand()*1.0/RAND_MAX;
	cout<<"vector["<<i<<"]="<<v4[i]<<endl;
	}
    // Вектор из элементов типа float
    // С неопределенным размером
    vector<float> v2;

    // Вектор, состоящий из 10 элементов типа int
    // По умолчанию все элементы заполняются нулями
    vector<int> v3(10, 0);

    return 0;
}
