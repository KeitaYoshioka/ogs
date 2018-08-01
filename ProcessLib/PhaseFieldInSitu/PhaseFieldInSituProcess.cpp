/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldInSituProcess.h"
#include "PhaseFieldInSituProcess-impl.h"

namespace ProcessLib
{
namespace PhaseFieldInSitu
{
template class PhaseFieldInSituProcess<2>;
template class PhaseFieldInSituProcess<3>;

}  // namespace PhaseFieldInSitu
}  // namespace ProcessLib
