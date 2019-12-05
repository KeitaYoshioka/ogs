/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DPHMPhaseFieldProcess.h"
#include "DPHMPhaseFieldProcess-impl.h"

namespace ProcessLib
{
namespace DPHMPhaseField
{
template class DPHMPhaseFieldProcess<2>;
template class DPHMPhaseFieldProcess<3>;

}  // namespace DPHMPhaseField
}  // namespace ProcessLib
