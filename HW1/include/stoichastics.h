#ifndef STOICHASTICS_H
#define STOICHASTICS_H

#include <vector>
#include <random>

namespace stoic{
    using namespace std;

    template<typename T>
    float mean(vector<T> v){
        return accumulate(v.begin(), v.end(), 0) / (float) v.size();
    };

    template<typename T>
    float var(vector<T> v){
        float v_ = mean(v);
        return inner_product(v.begin(), v.end(), v.begin(), 0.0) / v.size() - v_ * v_;
    }

    template<typename T>
    float sigma(vector<T> v){
        return sqrt(var(v));
    }

    vector<float> vec_r(int length, float low = 0.0f, float high = 1.0f){
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(low, high);
        vector<float> ret(length);
        for (int i=0; i<length; ++i) {
            ret.at(i) = distribution(generator);
        }
        return ret;
    }

    vector<float> vec_rn(int length, float mean = 0.0f, float dev = 1.0f){
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(mean, dev);
        vector<float> ret(length);
        for (int i=0; i<length; ++i) {
            ret.at(i) = distribution(generator);
        }
        return ret;
    }
};

#endif