/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroPhaseFieldProcess.h"
#include "HydroPhaseFieldProcess-impl.h"

namespace ProcessLib
{
namespace HydroPhaseField
{
template class HydroPhaseFieldProcess<2>;
template class HydroPhaseFieldProcess<3>;

}  // namespace HydroPhaseField
}  // namespace ProcessLib
