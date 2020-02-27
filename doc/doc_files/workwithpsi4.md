Working with Psi4 plugins     {#workwithpsi4}
==========================

[TOC]

@section pluginnote Short note on plugins

Psi4 plugins are 
[module libraries](https://cmake.org/cmake/help/v3.16/command/add_library.html?highlight=module%20library)
that are not intended to be linked to another target entity.
However, modules are created mainly because they can be dynamically
imported at runtime by consuming programs. The 
[input file](https://github.com/SinaMostafanejad/v2rdm_casscf/blob/master/tests/v2rdm_casscf_pdft/input.dat)
that we created to show the present example, is such a case where we import both
v2rdm\_casscf and mcpdft modules. At lines 51 and 53 of this input file, we call
two functions provided by these modules:

\code{.py}
en,wfn=energy('v2rdm-casscf',return_wfn=True)

energy('mcpdft',ref_wfn=wfn)
\endcode

The energy('v2rdm-casscf') call generates <i>wfn</i> object needed by the second
energy('mcpdft') call. Be sure to add the corresponding <i>return_wfn</i> and
<i>ref_wfn</i> arguments to these function calls, respectively as shown by the example.

@section psi4openrdm How OpenRDM connects with Psi4 plugins

There are at least two ways that you can choose to let OpenRDM talk to your
Psi4 plugin:

* Using python commands in your input files
* Using CMake commands in the CMakeList.txt file within plugin's source directory

Take a glance at a sample input file from our
[v2RDM-CASSCF plugin fork](https://github.com/SinaMostafanejad/v2rdm_casscf)
that imports the <i>MC-PDF</i> module from OpenRDM.

Please note that before adopting one of the following routes, make sure that
both Psi4 and OpenRDM are installed on your system in the aforementioned order.

@section usingpython Using python

Assuming that you have installed OpenRDM in <i>/home/OpenRDM</i> directory,
include the following command in the begining of your Psi4 input file and 
then import the MCPDFT module:

\code{.py}
sys.path.insert(0, '/home/OpenRDM/build/stage')
import mcpdft
\endcode

Check the [input file](https://github.com/SinaMostafanejad/v2rdm_casscf/blob/master/tests/v2rdm_casscf_pdft/input.dat)
in our [v2RDM-CASSCF plugin fork](https://github.com/SinaMostafanejad/v2rdm_casscf)
for more details.

@section usingcmake Using CMake

If the infrastructure of your software is managed by the CMake build system
generator, we recommend that your first take a glance at \ref howtouse section.
After finding out about what you need to add to the CMakeLists.txt file residing
in the main source directory of your project, CMake's internal
<i>find\_package()</i> takes care of finding and incorporating all necessary 
information about the OpenRDM by interrogating the openrdm::openrdm 
[IMPORTED target](https://cmake.org/cmake/help/v3.17/manual/cmake-buildsystem.7.html?highlight=imported#imported-targets).
