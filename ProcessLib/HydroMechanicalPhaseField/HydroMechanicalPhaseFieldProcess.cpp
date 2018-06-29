/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "HydroMechanicalPhaseFieldProcess.h"
#include "HydroMechanicalPhaseFieldProcess-impl.h"

namespace ProcessLib
{
namespace HydroMechanicalPhaseField
{
template class HydroMechanicalPhaseFieldProcess<2>;
template class HydroMechanicalPhaseFieldProcess<3>;

}  // namespace HydroMechanicalPhaseField
}  // namespace ProcessLib
