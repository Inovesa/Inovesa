#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#include "IO/GUI/GUIElement.hpp"
#include "IO/Display.hpp"

namespace vfps
{

class Plot3DColormap : public GUIElement
{
public:
	Plot3DColormap();

	~Plot3DColormap();

	void createTexture(PhaseSpace* mesh);

	void delTexture();

	void draw();

private:
    GLuint vertexbuffer;
    GLuint uvbuffer;

    GLuint UVArrayID;
    GLuint VertexArrayID;

    GLuint Texture;
    GLuint TextureID;
};

} // namespace vfps

#endif // PLOT2D_HPP
