#include "Edge.h"

using namespace std;

Edge::Edge(int v1, int v2, int f1, int f2) : v1(min(v1, v2)), v2(max(v1, v2)), f1(min(f1, f2)), f2(max(f1, f2)) {

}
