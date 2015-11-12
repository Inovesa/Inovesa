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

    void updateLine(PhaseSpace* mesh);

    void draw();

private:
    size_t _npoints;

    std::vector<float> _line;

    GLuint vertexbuffer;
    GLuint position;
};

} // namespace vfps

#endif // PLOT2DLINE_HPP
