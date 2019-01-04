/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldExternalProcess.h"
#include "PhaseFieldExternalProcess-impl.h"

namespace ProcessLib
{
namespace PhaseFieldExternal
{
template class PhaseFieldExternalProcess<2>;
template class PhaseFieldExternalProcess<3>;

}  // namespace PhaseFieldExternal
}  // namespace ProcessLib
