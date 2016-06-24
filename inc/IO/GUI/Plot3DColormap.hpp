/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#ifndef PLOT2D_HPP
#define PLOT2D_HPP

#ifdef INOVESA_USE_GUI

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

#endif // INOVESA_USE_GUI

#endif // PLOT2D_HPP
