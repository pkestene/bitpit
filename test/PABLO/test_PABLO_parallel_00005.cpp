/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include <mpi.h>

#include "bitpit_PABLO.hpp"

using namespace bitpit;

/*!
* Subtest 001
*
* Testing adaptive mesh refinement (AMR) with periodic boundary conditions.
*/
int subtest_001()
{
    /**<Instantation of a 3D pablo uniform object.*/
    PabloUniform pablo11(3);

    /** Set Periodic boundary conditions */
    pablo11.setPeriodic(0);
    pablo11.setPeriodic(2);
    pablo11.setPeriodic(4);

    /**<Refine globally one level and write the octree.*/
    pablo11.adaptGlobalRefine();

    /**<Refine globally one level and write the octree.*/
    pablo11.adaptGlobalRefine();

    /**<Partition the tree.*/
    pablo11.loadBalance();

    /**<Write the tree.*/
    pablo11.updateConnectivity();
    pablo11.write("pablo_parallel_00005");

    /**<Extract and print face neighbors.*/
    uint32_t nOct = pablo11.getNumOctants();
    std::cout << "Number of Octants : " << nOct << std::endl;
    std::cout << "Extracting the four face neighbours of each Octant " << std::endl;
    for (uint32_t iOct = 0; iOct < nOct; ++iOct) {
        std::vector<uint32_t> neigh;
        std::vector<bool> isGhost;
        std::cout << ". Octant index : " << iOct << std::endl;
        for (uint8_t iFace=0; iFace<pablo11.getNfaces(); ++iFace) {
            pablo11.findNeighbours(iOct, iFace, 1, neigh, isGhost);
            std::cout << " - For face " << (int)iFace << "; " << neigh.size() << " neighbours: [ ";
            for (auto iNeigh : neigh) {
	      std::cout << iNeigh << " " << (isGhost[0] ? "(ghost)" : "");
            }
            std::cout << "]" << std::endl;
        }
    }

    /**<Extract and print vertex neighbors .*/
    std::cout << "Extracting the four vertex neighbours of each Octant " << std::endl;
    for (uint32_t iOct = 0; iOct < nOct; ++iOct) {
        std::vector<uint32_t> neigh;
        std::vector<bool> isGhost;
        std::cout << ". Octant index : " << iOct << std::endl;
        for (uint8_t iVertex=0; iVertex<pablo11.getNnodes(); ++iVertex) {
            pablo11.findNeighbours(iOct, iVertex, 3, neigh, isGhost);
            std::cout << " - For vertex " << (int)iVertex << "; " << neigh.size() << " neighbours: [ ";
            for (auto iNeigh: neigh) {
                std::cout << iNeigh << " " << (isGhost[0] ? "(ghost)" : "");
            }
            std::cout << "]" << std::endl;
        }
    }

    /**<Extract and print edge neighbors .*/
    std::cout << "Extracting the four edge neighbours of each Octant " << std::endl;
    for (uint32_t iOct=0; iOct < nOct; ++iOct) {
        std::vector<uint32_t> neigh;
        std::vector<bool> isGhost;
        std::cout << ". Octant index : " << iOct << std::endl;
        for (uint8_t iEdge = 0; iEdge<pablo11.getNedges(); ++iEdge) {
            pablo11.findNeighbours(iOct, iEdge, 2, neigh, isGhost);
            std::cout << " - For edge " << (int)iEdge << "; " << neigh.size() << " neighbours: [ ";
            for (auto iNeigh: neigh) {
                std::cout << iNeigh << " " << (isGhost[0] ? "(ghost)" : "");
            }
            std::cout << "]" << std::endl;
        }
    }

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);

    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Run the subtests
    log::cout() << "Testing parallel AMR with periodic boundary conditions." << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    MPI_Finalize();
}
