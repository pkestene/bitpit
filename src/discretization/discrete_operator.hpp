/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#ifndef __BITPIT_DISCRETE_OPERATOR_HPP__
#define __BITPIT_DISCRETE_OPERATOR_HPP__

#include <vector>

#include "bitpit_LA.hpp"

#include "stencil.hpp"

namespace bitpit {

class DiscreteOperator : protected SystemSolver {

public:
#if BITPIT_ENABLE_MPI==1
    DiscreteOperator(MPI_Comm communicator, bool debug = false);
    DiscreteOperator(MPI_Comm communicator, long nUnknowns, long nNZ = 0, bool debug = false);
#else
    DiscreteOperator(bool debug = false);
    DiscreteOperator(long nUnknowns, long nMaximumNZ = 0, bool debug = false);
#endif

    void initialize(long nUnknowns, long nMaximumNZ);
    void clear(bool release = false);

    void appendStencil(const StencilScalar &stencil);

    void solve();
    void solve(const std::vector<double> &rhs, std::vector<double> *solution);

    using SystemSolver::getRHSRawPtr;
    using SystemSolver::getSolutionRawPtr;

protected:
    SparseMatrix m_matrix;
    std::vector<double> m_constants;

    void _initialize(long nUnknowns, long nMaximumNZ);

};

}

#endif
