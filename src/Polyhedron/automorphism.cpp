#include "Polyhedron.h"
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

struct Automorphism {
    std::vector<int> vertex_map;  // vertex_map[src] = dst
    std::vector<int> face_map;    // face_map[src] = dst
};

void Polyhedron::compute_automorphisms() {
    fix_winding();

    auto& f2v = get_sorted_f2v();
    auto& v2f = get_sorted_v2f();
    auto& v2v = get_v2v();

    const int V = v_count;
    const int F = f_count;

    // --- Precompute signatures for pruning ---
    // Vertex sig: sorted tuple of (face_size) for each incident face
    // Face sig: sorted tuple of (vertex_degree) for each incident vertex
    auto make_v_sig = [&](int v) -> std::vector<int> {
        std::vector<int> sig;
        sig.reserve(v2f[v].size());
        for (int f : v2f[v]) sig.push_back((int)f2v[f].size());
        std::sort(sig.begin(), sig.end());
        return sig;
    };

    auto make_f_sig = [&](int f) -> std::vector<int> {
        std::vector<int> sig;
        sig.reserve(f2v[f].size());
        for (int v : f2v[f]) sig.push_back((int)v2v[v].size());
        std::sort(sig.begin(), sig.end());
        return sig;
    };

    std::vector<std::vector<int>> v_sigs(V), f_sigs(F);
    for (int v = 0; v < V; v++) v_sigs[v] = make_v_sig(v);
    for (int f = 0; f < F; f++) f_sigs[f] = make_f_sig(f);

    // --- Precompute: for a directed edge (va, vb), find the face that has
    //     this edge in its winding (i.e., va then vb consecutively) ---
    // We'll need: given two vertices sharing an edge, find the two faces.
    // Given face f contains edge (va, vb) in winding order,
    // the OTHER face across that edge has (vb, va) in its winding.

    // Helper: find index of val in vec
    auto index_of = [](const std::vector<int>& vec, int val) -> int {
        for (int i = 0; i < (int)vec.size(); i++) {
            if (vec[i] == val) return i;
        }
        return -1;
    };

    // For face fs, edge at position i is (f2v[fs][i], f2v[fs][(i+1)%k]).
    // The adjacent face across that edge: find the face (other than fs)
    // that contains both vertices. We can find it via v2f intersection.
    // Cache this as face_adj[fs][i] = adjacent face index across edge i of fs.
    std::vector<std::vector<int>> face_adj(F);
    for (int f = 0; f < F; f++) {
        int k = (int)f2v[f].size();
        face_adj[f].resize(k);
        for (int i = 0; i < k; i++) {
            int va = f2v[f][i];
            int vb = f2v[f][(i + 1) % k];
            // Find face containing edge (vb, va) — the reversed winding
            int adj = -1;
            for (int ff : v2f[va]) {
                if (ff == f) continue;
                if (index_of(f2v[ff], vb) != -1) {
                    adj = ff;
                    break;
                }
            }
            face_adj[f][i] = adj;
        }
    }

    // --- Main search ---
    // Fix root flag: face 0, with v0 = f2v[0][0] at position 0
    const int f0 = 0;
    const int k0 = (int)f2v[f0].size();

    std::vector<Automorphism> results;

    // Try every target flag (ft, vt) where vt is at some position in f2v[ft]
    for (int ft = 0; ft < F; ft++) {
        // Quick prune: face sizes must match
        if ((int)f2v[ft].size() != k0) continue;
        // Signature prune
        if (f_sigs[ft] != f_sigs[f0]) continue;

        for (int vt_pos = 0; vt_pos < k0; vt_pos++) {
            int vt = f2v[ft][vt_pos];

            // Vertex signature prune
            if (v_sigs[vt] != v_sigs[f2v[f0][0]]) continue;

            // Attempt BFS extension
            std::vector<int> vmap(V, -1);
            std::vector<int> fmap(F, -1);
            bool valid = true;

            // Lambda to set vertex mapping with conflict check
            auto set_vmap = [&](int src, int dst) -> bool {
                if (vmap[src] == dst) return true;
                if (vmap[src] != -1) return false; // conflict
                // Also check injectivity: no other src maps to dst
                // For speed, we use a reverse map
                vmap[src] = dst;
                return true;
            };

            // Seed: align f0 → ft, with vertex rotation offset
            fmap[f0] = ft;
            for (int i = 0; i < k0; i++) {
                int vs = f2v[f0][i];
                int vd = f2v[ft][(vt_pos + i) % k0];
                if (!set_vmap(vs, vd)) { valid = false; break; }
            }
            if (!valid) continue;

            // BFS through face adjacency
            std::queue<int> queue;
            queue.push(f0);

            while (!queue.empty() && valid) {
                int fs = queue.front();
                queue.pop();

                int ks = (int)f2v[fs].size();
                int fd = fmap[fs];
                int kd = (int)f2v[fd].size();
                // ks == kd guaranteed if we got here

                for (int i = 0; i < ks; i++) {
                    int va_s = f2v[fs][i];
                    int vb_s = f2v[fs][(i + 1) % ks];

                    int va_d = vmap[va_s];
                    int vb_d = vmap[vb_s];
                    // These should both be set since we aligned this face

                    // Adjacent face in source across this edge
                    int f_adj_s = face_adj[fs][i];
                    if (f_adj_s == -1) { valid = false; break; }

                    // Adjacent face in destination across mapped edge
                    // Edge (va_d, vb_d) belongs to face fd in winding order.
                    // The adjacent face across it has (vb_d, va_d) in winding.
                    // Find that face:
                    int f_adj_d = -1;
                    for (int ff : v2f[va_d]) {
                        if (ff == fd) continue;
                        if (index_of(f2v[ff], vb_d) != -1) {
                            f_adj_d = ff;
                            break;
                        }
                    }
                    if (f_adj_d == -1) { valid = false; break; }

                    if (fmap[f_adj_s] != -1) {
                        // Already mapped — check consistency
                        if (fmap[f_adj_s] != f_adj_d) { valid = false; break; }
                        continue;
                    }

                    // Size check
                    if (f2v[f_adj_s].size() != f2v[f_adj_d].size()) {
                        valid = false; break;
                    }

                    fmap[f_adj_s] = f_adj_d;
                    queue.push(f_adj_s);

                    // Align vertices of f_adj_s → f_adj_d
                    // In f_adj_s, the shared edge appears as (..., vb_s, va_s, ...)
                    // because winding reverses across an edge.
                    // In f_adj_d, the shared edge appears as (..., vb_d, va_d, ...)
                    int k_adj = (int)f2v[f_adj_s].size();
                    int off_s = index_of(f2v[f_adj_s], vb_s);
                    int off_d = index_of(f2v[f_adj_d], vb_d);

                    if (off_s == -1 || off_d == -1) { valid = false; break; }

                    for (int j = 0; j < k_adj; j++) {
                        int vs = f2v[f_adj_s][(off_s + j) % k_adj];
                        int vd = f2v[f_adj_d][(off_d + j) % k_adj];
                        if (!set_vmap(vs, vd)) { valid = false; break; }
                    }
                }
            }

            if (!valid) continue;

            // Verify completeness — all vertices and faces mapped
            bool complete = true;
            for (int v = 0; v < V; v++) {
                if (vmap[v] == -1) { complete = false; break; }
            }
            if (!complete) continue;

            // Injectivity check on vertex map
            {
                std::vector<bool> used(V, false);
                bool injective = true;
                for (int v = 0; v < V; v++) {
                    if (used[vmap[v]]) { injective = false; break; }
                    used[vmap[v]] = true;
                }
                if (!injective) continue;
            }

            results.push_back({std::move(vmap), std::move(fmap)});
        }
    }
    //print results
    std::cout << "Found " << results.size() << " automorphisms:" << std::endl;
    for (const auto& auto_ : results) {
        std::cout << "Vertex map: ";
        for (int v = 0; v < V; v++) {
            std::cout << auto_.vertex_map[v] << " ";
        }
        std::cout << "\nFace map: ";
        for (int f = 0; f < F; f++) {
            std::cout << auto_.face_map[f] << " ";
        }
        std::cout << "\n---\n" << endl;
    }
}