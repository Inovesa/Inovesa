#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#include "IO/GUI/GUIElement.hpp"
#include "IO/Display.hpp"

namespace vfps
{

class Plot3DColormap : public GUIElement
{
public:
	Plot3DColormap(vfps::meshdata_t maxvalue=1);

	~Plot3DColormap();

	void createTexture(PhaseSpace* mesh);

	void delTexture();

	void draw();

private:
    GLuint vertexbuffer;
    GLuint uvbuffer;

    GLuint vertexUV;
    GLuint position;

    GLuint textureID;
    GLuint textureSampler;

    vfps::meshdata_t maxValue;
};

} // namespace vfps

#endif // PLOT2D_HPP
