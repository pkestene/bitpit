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

#include "discrete_operator.hpp"

namespace bitpit {

/*!
* \ingroup discretization
* \class DiscreteOperator
*
* The DiscreteOperator handles discrete operators.
*/

/*!
* Constuctor
*/
#if BITPIT_ENABLE_MPI==1
/*!
* \param communicator is the communicator that will be used
*/
#endif
/*!
* \param debug if this parameter is set to true, debug informations will be
* printed when solving the system
*/
#if BITPIT_ENABLE_MPI==1
DiscreteOperator::DiscreteOperator(MPI_Comm communicator, bool debug)
    : SystemSolver(communicator, debug), m_matrix(communicator)
#else
DiscreteOperator::DiscreteOperator(bool debug)
    : SystemSolver(debug), m_matrix()
#endif
{
}

/*!
* Constuctor
*/
#if BITPIT_ENABLE_MPI==1
/*!
* \param communicator is the communicator that will be used
*/
#endif
/*!
* \param nUnknowns is the number of unknowns
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
* \param debug if this parameter is set to true, debug informations will be
* printed when solving the system
*/
#if BITPIT_ENABLE_MPI==1
DiscreteOperator::DiscreteOperator(MPI_Comm communicator, long nUnknowns, long nNZ, bool debug)
    : SystemSolver(communicator, debug), m_matrix(communicator, nUnknowns, nUnknowns, nNZ)
#else
DiscreteOperator::DiscreteOperator(long nUnknowns, long nNZ, bool debug)
    : SystemSolver(debug), m_matrix(nUnknowns, nUnknowns, nNZ)
#endif
{
    _initialize(nUnknowns, nNZ);
}

/*!
* Initialize the discrete operator.
*
* \param nUnknowns is the number of unknowns
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
void DiscreteOperator::initialize(long nUnknowns, long nNZ)
{
    _initialize(nUnknowns, nNZ);
}

/*!
* Internal function to initialize the discrete operator.
*
* \param nUnknowns is the number of unknowns
* \param nNZ is the number of non-zero elements that the matrix will contain.
* This is just an optional hint. If the actual number of non-zero elements
* turns out to be greater than the provided value, the initialization of the
* matrix pattern will be slower because reallocation of internal data may be
* needed
*/
void DiscreteOperator::_initialize(long nUnknowns, long nNZ)
{
    m_matrix.initialize(nUnknowns, nUnknowns, nNZ);
    m_constants.resize(nUnknowns);
}

/*!
* Clear the discrete operator
*
* \param release if it's true the memory hold by the discrete operator will be
* released, otherwise the discrete operator will be cleared but its memory will
* not be relased
 */
void DiscreteOperator::clear(bool release)
{
    SystemSolver::clear();

    m_matrix.clear(release);
    if (release) {
        std::vector<double>().swap(m_constants);
    } else {
        m_constants.clear();
    }
}

/*!
* Append a stencil.
*
* \param stencil is the stencil that wil be appended
*/
void DiscreteOperator::appendStencil(const StencilScalar &stencil)
{
    m_matrix.addRow(stencil.getPattern().getItemCount(), stencil.getPattern().data(), stencil.getWeights().data());
    m_constants.push_back(stencil.getConstant());
}

/*!
* Solve the system
*/
void DiscreteOperator::solve()
{
    // Initialzie the system
    if (isInitialized()) {
        initialize(m_matrix.getRowCount(), m_matrix.getNZCount());
    }

    // Subtract constant terms to the RHS
    long nUnknowns = m_matrix.getRowCount();
    double *raw_rhs = getRHSRawPtr();
    for (long i = 0; i < nUnknowns; ++i) {
        raw_rhs[i] -= m_constants[i];
    }
    restoreRHSRawPtr(raw_rhs);

    // Solve the system
    SystemSolver::solve();
}

/*!
* Solve the system
*
* \param rhs is the right-hand-side of the system
* \param solution in input should contain the initial solution, on output it
* contains the solution of the linear system
*/
void DiscreteOperator::solve(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Fills the vectors
    vectorsFill(rhs, solution);

    // Solve the system
    solve();

    // Export the solution
    vectorsExport(solution);
}

}
