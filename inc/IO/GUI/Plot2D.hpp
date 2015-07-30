#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#include "IO/GUI/GUIElement.hpp"
#include "IO/Display.hpp"

namespace vfps
{

class Plot2D : public GUIElement
{
public:
	Plot2D();

	~Plot2D();

	void createTexture(PhaseSpace* mesh);

	void delTexture();

	void draw();
};

} // namespace vfps

#endif // PLOT2D_HPP
