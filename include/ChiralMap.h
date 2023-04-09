#ifndef CHIRALMAP_H
#define CHIRALMAP_H
#include <tuple>
#include <map>
#include <unordered_map>
#include <vector>
#include <functional>

using namespace std;

//this class puts 2 vertices on each edge, and

struct ChiralKey
{
    ChiralKey(int v1_, int v2_, int v3_){
        assign(v1_,v2_,v3_);
    }

    void assign(int v1_, int v2_, int v3_){
        if(v1_ == v3_){
            if(v2_ > v1_){
                v1 = v1_;
                v2 = v3_;
                v3 = v2_;
            }else{
                v1 = v2_;
                v2 = v3_;
                v3 = v1_;
            }
        }else{ //v2 == v3
            if(v2_ > v1_){
                v1 = v1_;
                v2 = v2_;
                v3 = v3_;
            }else{
                v1 = v2_;
                v2 = v3_;
                v3 = v1_;
            }
        }
    }
    int v1, v2, v3;
};

namespace std
{
    template<>
    struct hash<ChiralKey>
    {
        std::size_t operator()(const ChiralKey& key) const
        {
            //based on jenkins hash function.
            //the magic number is the reciprocal of the golden ratio
            //bob jenkins:
            //"The golden ratio really is an arbitrary value. Its purpose is to avoid mapping all zeros to all zeros."
            //but tbh we don't have that problem anyways.
            size_t seed = std::hash<int>()(key.v1);
            seed ^= hash<int>()(key.v2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= hash<int>()(key.v3) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            return seed;
        }
    };

    template<>
    struct equal_to<ChiralKey>
    {
        bool operator()(const ChiralKey& lhs, const ChiralKey& rhs) const
        {
            return (lhs.v1 == rhs.v1) && (lhs.v2 == rhs.v2) && (lhs.v3 == rhs.v3);
        }
    };
}

typedef unordered_map<ChiralKey,int> ChiralMap;

#endif
