// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
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
