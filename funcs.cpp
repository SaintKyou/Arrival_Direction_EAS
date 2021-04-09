#include "funcs.h"
#include <vector>
#include <algorithm>

int dec2bit(short num)
{
    int bin = 0, k =1;
    while(num){
    bin+=(num%2)*k;
    k*=10;
    num/=2;
}
    return bin;
}
int get_dig(int num)
{
    int count{}, k = 1, i = 0;
    while (i<9) {
       if(num%10==1){
        ++count;
        }
       ++i;
       num/=10;
    }
    return count;
}
int get_num(int num)
{
    int count{}, k = 1, i = 0;
    while (i<9) {
       if(num!=0){
        ++count;
        }
       ++i;
       num/=10;
    }
    return count;
}
std::vector<unsigned int> make_v(int num){
    std::vector<unsigned int> v;
    while(num/10!=0)
    {
        v.push_back(num%10);
        num/=10;
    }
    v.push_back(1);
    while(v.size()<9){
        v.push_back(0);
    }
    reverse(v.begin(),v.end());
    return v;
}

std::vector<unsigned int> make_v2(int num){
    std::vector<unsigned int> v;
    while(num/10!=0)
    {
        v.push_back(num%10);
        num/=10;
    }
    v.push_back(1);
    while(v.size()<4){
        v.push_back(0);
    }
    reverse(v.begin(),v.end());
    return v;
}

std::vector<long> norm_vec (std::vector<long> &a)
{
    long min = a[0];
    for(int i = 1; i < a.size(); ++i){
        if(a[i] < min) min = a[i];
    }
    for(int i = 0; i < a.size(); ++i) a[i]-=min;

}
