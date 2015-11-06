#ifndef PLOT2DLINE_HPP
#define PLOT2DLINE_HPP

#include <vector>

#include "IO/GUI/GUIElement.hpp"

namespace vfps
{

class Plot2DLine : public GUIElement
{
public:
    Plot2DLine();

    ~Plot2DLine();

    void createLine(PhaseSpace* mesh);

    void delLine();

    void draw();

private:
    size_t _npoints;

    std::vector<float> _line;
};

} // namespace vfps

#endif // PLOT2DLINE_HPP
