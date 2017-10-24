/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include "ProcessLib/Process.h"

namespace BaseLib
{
class ConfigTree;
}
namespace MeshLib
{
class Mesh;
}
namespace ProcessLib
{
class AbstractJacobianAssembler;
struct ParameterBase;
class Process;
class ProcessVariable;
}

namespace ProcessLib
{
namespace PhaseFieldStaggered
{
std::unique_ptr<Process> createPhaseFieldStaggeredProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace PhaseFieldStaggered
}  // namespace ProcessLib
