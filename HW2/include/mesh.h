#pragma once

#include <string>
#include <vector>

namespace euler {

    constexpr int NG = 2; // number of ghost layers on every side of the block

    /*
       Single-block structured O-grid.

       Nodes:   physical i = 0..ni-1 (around the airfoil, TE -> lower -> LE -> upper -> TE)
                physical j = 0..nj-1 (j = 0 wall, j = nj-1 farfield)
                node (0, j) and node (ni-1, j) are the same point (cut at the TE).
       Cells:   cell (i, j) has corners (i,j), (i+1,j), (i+1,j+1), (i,j+1);
                physical i = 0..nic-1, j = 0..njc-1.
       Ghosts:  every index extends by NG on both sides, e.g. i = -NG..ni-1+NG,
                so neighbours of boundary cells are reached with plain i-1, j-1.
       Storage: one flat std::vector per quantity, i contiguous,
                so inner loops over i walk memory sequentially.
    */
    struct Mesh {
        int ni = 0;      // physical nodes in i 
        int nj = 0;      // physical nodes in j (wall to farfield)
        int nic = 0;     // physical cells in i (= ni - 1)
        int njc = 0;     // physical cells in j (= nj - 1)
        int ni_tot = 0;  // stored nodes in i, ghosts included (= ni + 2*NG)
        int nj_tot = 0;  // stored nodes in j, ghosts included (= nj + 2*NG)
        int nic_tot = 0; // stored cells in i, ghosts included (= nic + 2*NG)
        int njc_tot = 0; // stored cells in j, ghosts included (= njc + 2*NG)

        std::vector<double> x; // node x coordinates, size ni_tot * nj_tot (ghosts = NaN until filled)
        std::vector<double> y; // node y coordinates, size ni_tot * nj_tot (ghosts = NaN until filled)

        /*
           Metrics (Euler2D slide 21), filled by compute_metrics().
           A face is stored at the index of its first node, node(i, j):
             i-face (i, j): edge node(i,j) -> node(i,j+1),  i = 0..ni-1, j = 0..njc-1,
                            "left" cell (i-1, j), "right" cell (i, j)
             j-face (i, j): edge node(i,j) -> node(i+1,j),  i = 0..nic-1, j = 0..nj-1,
                            "left" cell (i, j-1), "right" cell (i, j)
           Normals are unit vectors from left to right (towards +i or +j).
           At the boundaries the left or right cell is a ghost cell.
           Unused slots (ghosts) stay NaN.
        */
        std::vector<double> vol; // cell area Omega (the 2D "volume"), size nic_tot * njc_tot
        std::vector<double> nix; // i-face unit normal, x component, size ni_tot * nj_tot
        std::vector<double> niy; // i-face unit normal, y component, size ni_tot * nj_tot
        std::vector<double> si;  // i-face length Delta S, size ni_tot * nj_tot
        std::vector<double> njx; // j-face unit normal, x component, size ni_tot * nj_tot
        std::vector<double> njy; // j-face unit normal, y component, size ni_tot * nj_tot
        std::vector<double> sj;  // j-face length Delta S, size ni_tot * nj_tot

        // Flat index of node (i, j); i, j may be negative (ghosts).
        int node(int i, int j) const { return (j + NG) * ni_tot + (i + NG); }

        // Flat index of cell (i, j); same convention. Used from step 2 on.
        int cell(int i, int j) const { return (j + NG) * nic_tot + (i + NG); }
    };

    /*
       Read a single-block 2D plot3D file (Euler2D slide 18):
         nblock
         ni nj
         all x (i fastest, then j), then all y in the same order
    */
    Mesh read_plot3d(const std::string& filename);

    /*
       Compute the cell areas, face unit normals and face lengths (Euler2D slide 21)
       of the physical cells and faces. Throws if a cell area is not positive
       (cell corners not counter-clockwise, or a folded cell).
    */
    void compute_metrics(Mesh& m);

}
