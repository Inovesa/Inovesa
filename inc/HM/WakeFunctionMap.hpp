#ifndef WAKEFUNCTIONMAP_HPP
#define WAKEFUNCTIONMAP_HPP

#include "HM/WakeKickMap.hpp"

namespace vfps
{

class WakeFunctionMap : public WakeKickMap
{
public:
    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param wake wakefunction from -xsize to xsize-1, normalized in a way
     *             that plain multiplication with density gives an offset
     *
     * @todo interpolation when wake does not match PhaseSpace
     */
    WakeFunctionMap(PhaseSpace* in, PhaseSpace* out,
                    const meshindex_t xsize, const meshindex_t ysize,
                    const std::vector<std::pair<meshaxis_t,double>> wake,
                    const InterpolationType it);

    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param csrimpedance
     * @param it
     */
    WakeFunctionMap(PhaseSpace* in, PhaseSpace* out,
                    const meshindex_t xsize, const meshindex_t ysize,
                    const vfps::ElectricField* csr,
                    const InterpolationType it);

    ~WakeFunctionMap();

public:
    /**
     * @brief overloads KickMap::apply() to have a variable HeritageMap
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void apply();

    inline const meshaxis_t* getWakeFunction() const
        { return _wakefunction; }

private:


private:
    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param it
     */
    WakeFunctionMap(PhaseSpace* in, PhaseSpace* out,
                    const meshindex_t xsize, const meshindex_t ysize,
                    const InterpolationType it);

    /**
     * @brief _wakefunktion (normalized single particle) wake,
     *          sampled at 2*xsize positions [-xsize:+xsize]
     */
    meshaxis_t* const _wakefunction;

    const size_t _wakesize;
};

} // namespace VFPS

#endif // WAKEFUNCTIONMAP_HPP
