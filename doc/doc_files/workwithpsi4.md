Working with Psi4 plugins     {#workwithpsi4}
==========================

[TOC]

@section pluginnote A brief note on Psi4 plugins

Psi4 plugins are 
(shared) [module libraries](https://cmake.org/cmake/help/v3.16/command/add_library.html?highlight=module)
that are <b>not</b> intended to be linked to another target entity.
However, modules are created mainly because they can dynamically be
imported at runtime by the consuming programs. In our sample  
[input file](https://github.com/SinaMostafanejad/v2rdm_casscf/blob/master/tests/v2rdm_casscf_pdft/input.dat)
, we show how to import both <i>v2rdm\_casscf</i> and <i>mcpdft</i> modules.
This input file assumes that you have chosen the second way below (recommended) for 
connecting OpenRDM with your plugin. At lines 51 and 53 of this input file, we call
two functions that are provided by these modules:

\code{.py}
en,wfn=energy('v2rdm-casscf',return_wfn=True)

energy('mcpdft',ref_wfn=wfn)
\endcode

The <i>energy('v2rdm-casscf')</i> call generates a <i>wfn</i> object that is 
needed by the second <i>energy('mcpdft')</i> call. Be sure to add the corresponding
<i>return_wfn</i> and <i>ref_wfn</i> arguments to these function calls, respectively as shown.

@section psi4openrdm How does OpenRDM connect with Psi4 plugins?

There are at least two ways that you can adopt to let <b>OpenRDM</b> talk to your
Psi4 plugin:

* Using python commands in your input files
* Using CMake commands in the <i>CMakeList.txt</i> file within plugin's source directory

Before adopting one of these routes, make sure that
both Psi4 and <b>OpenRDM</b> are installed on your system in the aforementioned order
as the latter becomes dependent on the former in the present case.

When installing <b>OpenRDM</b>, turn on the CMake option <i>WITH_PSI4=ON</i>
(see [configure](https://github.com/SinaMostafanejad/OpenRDM/blob/master/configure)
file in <b>OpenRDM</b>'s main source directory for more details). Activating this
option requires that the CMake <i>psi4_DIR</i> variable in the configure script to be set
to the directory that contains the <i>psi4Config.cmake</i> file so that CMake's <i>find\_packge()</i>
function within <b>OpenRDM</b> is able to find the <i>psi4::core</i> IMPORTED target
and links it to <b>OpenRDM</b>.

@subsection usingpython Using python

Assuming that you have installed <b>OpenRDM</b> in <i>/home/OpenRDM</i> directory,
include the following command in the begining of your Psi4 input file and 
then import the MCPDFT module:

\code{.py}
sys.path.insert(0, '/home/OpenRDM/build/stage')
import mcpdft
\endcode

In the present example, the <i>mcpdft.so</i> module lives in
the <i>/home/<b>OpenRDM</b>/build/stage/lib</i> directory.

Check the [input file](https://github.com/SinaMostafanejad/v2rdm_casscf/blob/master/tests/v2rdm_casscf_pdft/input.dat)
in our [v2RDM-CASSCF plugin fork](https://github.com/SinaMostafanejad/v2rdm_casscf)
for more details.

@subsection usingcmake Using CMake

If the infrastructure of your software is managed by the CMake build system
generator, we recommend that your first take a glance at \ref howtouse section.
After finding out about what you need to add to the CMakeLists.txt file residing
in the main source directory of your project, CMake's internal
<i>find\_package()</i> takes care of finding and incorporating all necessary 
information about the <b>OpenRDM</b> by interrogating the <i>openrdm::openrdm</i> 
[IMPORTED target](https://cmake.org/cmake/help/v3.17/manual/cmake-buildsystem.7.html?highlight=imported#imported-targets).
