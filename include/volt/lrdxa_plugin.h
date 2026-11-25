#pragma once

#include <volt/plugin/plugin_main.h>
#include <volt/structures/crystal_structure_types.h>

#include <string>

namespace Volt::LrdxaPlugin {

Plugin::PluginDescriptor buildDescriptor();

// The crystal structure an option names; FCC, with a warning, for a name it does not know.
LatticeStructureType parseCrystalStructure(const std::string& name);

}
