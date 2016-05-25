#ifndef WAKEIMPEDANCEMAP_HPP
#define WAKEIMPEDANCEMAP_HPP

#include "SM/WakeKickMap.hpp"
#include "ElectricField.hpp"
#include "Impedance.hpp"

namespace vfps
{

class WakePotentialMap : public WakeKickMap
{
public:
    WakePotentialMap(PhaseSpace* in, PhaseSpace* out,
                     const meshindex_t xsize,
                     const meshindex_t ysize,
                     ElectricField* field,
                     const InterpolationType it,
                     bool interpol_clamp);

public:
    /**
     * @brief update implements WakeKickMap
     */
    void update();

private:
    ElectricField* _field;
};

} // namespace vfps

#endif // WAKEIMPEDANCEMAP_HPP
