// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#if INOVESA_USE_PNG == 1
#include <png++/png.hpp>
#endif
#include <iterator>

#include "defines.hpp"
#include "IO/Display.hpp"
#include "IO/HDF5File.hpp"
#include "PS/PhaseSpaceFactory.hpp"

#if INOVESA_USE_HDF5 == 1
std::unique_ptr<vfps::PhaseSpace>
vfps::makePSFromHDF5(std::string fname, int64_t startdiststep
                    , vfps::meshaxis_t qmin, vfps::meshaxis_t qmax
                    , vfps::meshaxis_t pmin, vfps::meshaxis_t pmax
                    , oclhptr_t oclh
                    , const double bunch_charge
                    , const double bunch_current
                    , double xscale, double yscale
                    )
{
    try {
        auto ps = HDF5File::readPhaseSpace(fname,qmin,qmax,pmin,pmax
                                          , oclh
                                          , bunch_charge,bunch_current,
                                           xscale,yscale,startdiststep);
        #if INOVESA_USE_OPENCL == 1
        ps->syncCLMem(OCLH::clCopyDirection::cpu2dev);
        #endif // INOVESA_USE_OPENCL
        return ps;
    } catch (const std::exception& ex) {
        std::cerr << "Error loading initial distribution from \""
                  << fname << "\":"
                  << ex.what() << std::endl;
    } catch (const H5::Exception& ex) {
        #if H5_VERS_MAJOR == 1 and H5_VERS_MINOR < 10
        ex.printError();
        #else
        ex.printErrorStack();
        #endif
        H5::Exception::clearErrorStack();
    } catch (...) {
        std::cerr << "Error loading initial distribution from \""
                  << fname << "\".";
    }
    return nullptr;
}
#endif // INOVESA_USE_HDF5

std::unique_ptr<vfps::PhaseSpace>
vfps::makePSFromPNG( std::string fname
                   , meshaxis_t qmin, meshaxis_t qmax
                   , meshaxis_t pmin, meshaxis_t pmax
                   , oclhptr_t oclh
                   , const double bunch_charge
                   , const double bunch_current
                   , double qscale, double pscale
                   )
{
    #if INOVESA_USE_PNG == 1
    // load pattern to start with
    png::image<png::gray_pixel_16> image;
    try {
        image.read(fname);
    } catch ( const png::std_error &e ) {
        std::cerr << e.what() << std::endl;
    } catch ( const png::error &e ) {
        std::cerr << "Problem loading " << fname
                  << ": " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Error loading initial distribution from \""
                  << fname << "\".";
    }

    if (image.get_width() == image.get_height()) {
        meshindex_t ps_size = image.get_height();

        std::vector<meshdata_t> data(ps_size*ps_size);

        for (unsigned int x=0; x<ps_size; x++) {
            for (unsigned int y=0; y<ps_size; y++) {
                data[x*ps_size+y] = image[ps_size-y-1][x]/float(UINT16_MAX);
            }
        }
        auto ps = std::make_unique<PhaseSpace>( ps_size
                                              , qmin, qmax, qscale
                                              , pmin, pmax, pscale
                                              , oclh
                                              , bunch_charge,bunch_current
                                              , 1U, 1
                                              , data.data());
        // normalize integral to 1
        ps->updateXProjection();
        ps->normalize();

        #if INOVESA_USE_OPENCL == 1
        ps->syncCLMem(OCLH::clCopyDirection::cpu2dev);
        #endif // INOVESA_USE_OPENCL
        std::stringstream imgsize;
        imgsize << ps_size;
        Display::printText("Read phase space (a="+imgsize.str()+" px).");
        return ps;
    } else {
        std::cerr << "Phase space has to be quadratic. Please adjust "
                  << fname << std::endl;
    }
#endif // INOVESA_USE_PNG
    return nullptr;
}

std::unique_ptr<vfps::PhaseSpace>
vfps::makePSFromTXT(std::string fname, int64_t ps_size
                   , vfps::meshaxis_t qmin, vfps::meshaxis_t qmax
                   , vfps::meshaxis_t pmin, vfps::meshaxis_t pmax
                   , oclhptr_t oclh
                   , const double bunch_charge, const double bunch_current
                   , double qscale, double pscale)
{
    auto ps = std::make_unique<PhaseSpace>( ps_size
                                          , qmin, qmax, qscale
                                          , pmin, pmax, pscale
                                          , oclh
                                          , bunch_charge,bunch_current
                                          , 1U
                                          );
    std::ifstream ifs;
    try {
        ifs.open(fname);
    } catch (const std::exception& ex) {
        std::cerr << "Error loading initial distribution from \""
                  << fname << "\":"
                  << ex.what() << std::endl;
    } catch (...) {
        std::cerr << "Error loading initial distribution from \""
                  << fname << "\".";
    }

    ifs.unsetf(std::ios_base::skipws);

    // count the newlines with an algorithm specialized for counting:
    size_t line_count = std::count(
        std::istream_iterator<char>(ifs),
        std::istream_iterator<char>(),
        '\n');

    ifs.setf(std::ios_base::skipws);
    ifs.clear();
    ifs.seekg(0,ifs.beg);

    while (ifs.good()) {
        float xf,yf;
        ifs >> xf >> yf;
        meshindex_t x = std::lround((xf/qmax+0.5f)*ps_size);
        meshindex_t y = std::lround((yf/pmax+0.5f)*ps_size);
        if (x < ps_size && y < ps_size) {
            (*ps)[0][x][y] += 1.0/line_count;
        }
    }
    ifs.close();

    // normalize integral to 1
    ps->normalize();
    #if INOVESA_USE_OPENCL == 1
    ps->syncCLMem(OCLH::clCopyDirection::cpu2dev);
    #endif // INOVESA_USE_OPENCL
    return ps;
}
