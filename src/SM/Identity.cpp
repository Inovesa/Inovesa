// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "SM/Identity.hpp"

vfps::Identity::~Identity() noexcept
#if INOVESA_ENABLE_CLPROFILING == 1
{
    saveTimings("IDM");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING
