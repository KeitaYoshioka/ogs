/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace PhaseFieldStaggered
{
struct PhaseFieldStaggeredProcessData
{
    PhaseFieldStaggeredProcessData(Parameter<double> const& residual_stiffness_,
                                   Parameter<double> const& crack_resistance_,
                                   Parameter<double> const& crack_length_scale_)
        : residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_)
    {
    }

    PhaseFieldStaggeredProcessData(PhaseFieldStaggeredProcessData&& other)
        : residual_stiffness(other.residual_stiffness),
          crack_resistance(other.crack_resistance),
          crack_length_scale(other.crack_length_scale),
          dt{other.dt},
          t{other.t}
    {
    }

    //! Copies are forbidden.
    PhaseFieldStaggeredProcessData(PhaseFieldStaggeredProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldStaggeredProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldStaggeredProcessData&&) = delete;

    Parameter<double> const& residual_stiffness;
    Parameter<double> const& crack_resistance;
    Parameter<double> const& crack_length_scale;
    double dt = 0;
    double t =0;
};

}  // namespace PhaseFieldStaggered
}  // namespace ProcessLib
