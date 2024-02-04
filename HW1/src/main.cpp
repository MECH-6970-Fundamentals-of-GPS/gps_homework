#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>

#include <matplot/matplot.h>
using namespace std;
using namespace matplot;

template<typename T>
vector<T> circshift(vector<T> vec, int index){
    if(index > 0){
        for(int i = 0; i < index; i++){
            vec.push_back(vec.front());
            vec.erase(vec.begin());
        }
    }else{
        for(int i = index; i < 0; i++){
            vec.insert(vec.begin(), vec.back());
            vec.pop_back();
        }
    }
    return vec;
}

template <typename T>
vector<T> corr(vector<T> v1){
    int N = v1.size();
    vector<T> ret;
    for(int i = 0; i < 2*N; i++){
        vector<T> shift = circshift(v1, i - N);
        T sum = 0;
        for(int j = 0; j < N; j++){
            sum = sum + (v1.at(j)*shift.at(j));
        }
        ret.insert(ret.begin(), sum / N);
    }
    return ret;
}

template <typename T>
vector<T> corr(vector<T> v1, vector<T> v2){
    int N = v1.size();
    vector<T> ret;
    for(int i = 0; i < 2*N; i++){
        vector<T> shift = circshift(v2, i - N);
        T sum = 0;
        for(int j = 0; j < N; j++){
            sum = sum + (v1.at(j)*shift.at(j));
        }
        ret.insert(ret.begin(), sum / N);
    }
    return ret;
}

int main() {
    int N = 1000;
    vector<float> r1(N);
    vector<float> r2(N);
    vector<int> idx(2*N);
    iota (begin(idx), end(idx), 0);
    for(int i = 0; i < r1.size(); i++){
        r1.at(i) = 2*ceil((rand() % 2) - 0.5) - 1;
        r2.at(i) = 2*ceil((rand() % 2) - 0.5) - 1;
    }
    vector<float> ac1 = corr(r1);
    vector<float> ac2 = corr(r2);
    vector<float> xc1 = corr(r1, r2);
    vector<float> xc2 = corr(r2, r1);

    figure();
    tiledlayout(2,1);
    nexttile();
    title("Autocorrelation of Random Sequence 1");
    xlabel("Index");
    ylabel("Correlation");
    plot(idx, ac1);

    nexttile();
    title("Autocorrelation of Random Sequence 2");
    xlabel("Index");
    ylabel("Correlation");
    plot(idx, ac2);

    figure();
    tiledlayout(2,1);
    nexttile();
    title("Cross Correlation of Random Sequence 1 w/ 2");
    xlabel("Index");
    ylabel("Correlation");
    plot(idx, xc1);

    nexttile();
    title("Cross Correlation of Random Sequence 2 w/ 1");
    xlabel("Index");
    ylabel("Correlation");
    plot(idx, xc2);

    show();
    return 0;
}