#ifndef EDGE_H
#define EDGE_H

#include <unordered_set>

using namespace std;

struct Edge
{
    int v1, v2; // vertex indices
    int f1, f2; // face indices

    Edge(int v1, int v2, int f1, int f2);
};

namespace std
{
    template<>
    struct hash<Edge>
    {
        std::size_t operator()(const Edge& edge) const
        {
            //based on jenkins hash function.
            //the magic number is the reciprocal of the golden ratio
            //bob jenkins:
            //"The golden ratio really is an arbitrary value. Its purpose is to avoid mapping all zeros to all zeros."
            //but tbh we don't have that problem anyways.
            size_t seed = std::hash<int>()(edge.v1);
            seed ^= std::hash<int>()(edge.v2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            //size_t seed = ((v1 + v2) * (v1 + v2 + 1)) / 2 + v2; //based on cantor pairing function. idk if better. nah.
            //i ran a bunch of tests of the cantor pairing function on finite rings up to size 32.
            //seems to be fairly random. but this doesnt say anything. using tried methods like jenkins' hash is more
            //reliable.
            return seed;
        }
    };

    template<>
    struct equal_to<Edge>
    {
        bool operator()(const Edge& lhs, const Edge& rhs) const
        {
            return (lhs.v1 == rhs.v1) && (lhs.v2 == rhs.v2);
        }
    };
}

typedef unordered_set<Edge> Edges;

#endif // EDGE_H
