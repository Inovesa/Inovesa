// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once


#include "SM/WakeKickMap.hpp"
#include "PS/ElectricField.hpp"
#include "Z/Impedance.hpp"

namespace vfps
{

class WakePotentialMap : public WakeKickMap
{
public:
    WakePotentialMap(std::shared_ptr<PhaseSpace> in
                    , std::shared_ptr<PhaseSpace> out
                    , ElectricField* field
                    , const InterpolationType it
                    , bool interpol_clamp
                    , oclhptr_t oclh
                    );

    ~WakePotentialMap() noexcept override;

public:
    /**
     * @brief update implements WakeKickMap
     */
    void update() override;

private:
    ElectricField* _field;
};

} // namespace vfps

