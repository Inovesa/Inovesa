#ifndef WAKEIMPEDANCEMAP_HPP
#define WAKEIMPEDANCEMAP_HPP

#include "HM/WakeKickMap.hpp"
#include "ElectricField.hpp"
#include "Impedance.hpp"

namespace vfps
{

class WakePotentialMap : public WakeKickMap
{
public:
    WakePotentialMap(PhaseSpace* in, PhaseSpace* out,
                     const meshindex_t xsize, const meshindex_t ysize,
                     ElectricField *field, const InterpolationType it);

public:
    /**
     * @brief overloads WakeKickMap::apply()
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void apply();

private:
    ElectricField* _field;
};

} // namespace vfps

#endif // WAKEIMPEDANCEMAP_HPP