#include <iostream>
#include "posit/posit"
#include "jack_settings.hpp"


#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>

#include "utils/outstream.hpp"


#include "Box.hpp"
#include "utils/Parameters.hpp"
#include "utils/utils.hpp"
#include "utils/mytimer.cpp"
#include "YAML_DOC.cpp"
#include "box_utils.hpp"
#include "driver.hpp"
#include "BoxPartition.hpp"


//The following macros should be specified as compile-macros in the
//makefile. They are defaulted here just in case...
#ifndef MINIFE_SCALAR
#define MINIFE_SCALAR positX;
#endif
#ifndef MINIFE_LOCAL_ORDINAL
#define MINIFE_LOCAL_ORDINAL int
#endif
#ifndef MINIFE_GLOBAL_ORDINAL
#define MINIFE_GLOBAL_ORDINAL int
#endif


// ************************************************************************

void add_params_to_yaml(YAML_Doc& doc, miniFE::Parameters& params);
void add_configuration_to_yaml(YAML_Doc& doc, int numprocs, int numthreads);
void add_timestring_to_yaml(YAML_Doc& doc);

//
//We will create a 'box' of size nx X ny X nz, partition it among processors,
//then call miniFE::driver which will use the partitioned box as the domain
//from which to assemble finite-element matrices into a global matrix and
//vector, then solve the linear-system using Conjugate Gradients.
//

int main(int argc, char** argv) {
    miniFE::Parameters params;
    miniFE::get_parameters(argc, argv, params);

    int numprocs = 1, myproc = 0;
    miniFE::initialize_mpi(argc, argv, numprocs, myproc);

    miniFE::timer_type start_time = miniFE::mytimer();

#ifdef MINIFE_DEBUG
    outstream(numprocs, myproc);
#endif

    //make sure each processor has the same parameters:
    miniFE::broadcast_parameters(params);


    Box global_box = { 0, params.nx, 0, params.ny, 0, params.nz };
    std::vector<Box> local_boxes(numprocs);

    box_partition(0, numprocs, 2, global_box, &local_boxes[0]);

    Box& my_box = local_boxes[myproc];

    MINIFE_GLOBAL_ORDINAL num_my_ids = miniFE::get_num_ids<MINIFE_GLOBAL_ORDINAL>(my_box);
    MINIFE_GLOBAL_ORDINAL min_ids = num_my_ids;

#ifdef HAVE_MPI
    MPI_Datatype mpi_dtype = miniFE::TypeTraits<MINIFE_GLOBAL_ORDINAL>::mpi_type();
  MPI_Allreduce(&num_my_ids, &min_ids, 1, mpi_dtype, MPI_MIN, MPI_COMM_WORLD);
#endif

    if (min_ids == 0) {
        std::cout<<"One or more processors have 0 equations. Not currently supported. Exiting."<<std::endl;

        miniFE::finalize_mpi();

        return 1;
    }

    std::ostringstream osstr;
    osstr << "miniFE." << params.nx << "x" << params.ny << "x" << params.nz;
#ifdef HAVE_MPI
    osstr << ".P"<<numprocs;
#endif
    osstr << ".";
    if (params.name != "") osstr << params.name << ".";

    YAML_Doc doc("miniFE", "Jack.0", ".", osstr.str());
    if (myproc == 0) {
        add_params_to_yaml(doc, params);
        add_configuration_to_yaml(doc, numprocs, params.numthreads);
        add_timestring_to_yaml(doc);
    }

    //Most of the program is performed in the 'driver' function, which is
    //templated on < Scalar, LocalOrdinal, GlobalOrdinal >.
    //To run miniFE with float instead of double, or 'long long' instead of int,
    //etc., change these template-parameters by changing the macro definitions in
    //the makefile or on the make command-line.

    int return_code = miniFE::driver< MINIFE_SCALAR, MINIFE_LOCAL_ORDINAL, MINIFE_GLOBAL_ORDINAL>(global_box, my_box, params, doc);

    miniFE::timer_type total_time = miniFE::mytimer() - start_time;

#ifdef MINIFE_REPORT_RUSAGE
    struct rusage get_mem;
   getrusage(RUSAGE_SELF, &get_mem);

   long long int rank_rss = get_mem.ru_maxrss;
   long long int global_rss = 0;
   long long int max_rss = 0;

#ifdef HAVE_MPI
   MPI_Reduce(&rank_rss, &global_rss, 1,
       	MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&rank_rss, &max_rss, 1,
       	MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
   if (myproc == 0) {
        doc.add("Global All-RSS (kB)", global_rss);
       	doc.add("Global Max-RSS (kB)", max_rss);
   }
#else
   doc.add("RSS (kB)", rank_rss);
#endif
#endif

    if (myproc == 0) {
        doc.add("Total Program Time",total_time);
        doc.generateYAML();
    }

    miniFE::finalize_mpi();

    return return_code;
}

void add_params_to_yaml(YAML_Doc& doc, miniFE::Parameters& params)
{
    doc.add("Global Run Parameters","");
    doc.get("Global Run Parameters")->add("dimensions","");
    doc.get("Global Run Parameters")->get("dimensions")->add("nx",params.nx);
    doc.get("Global Run Parameters")->get("dimensions")->add("ny",params.ny);
    doc.get("Global Run Parameters")->get("dimensions")->add("nz",params.nz);
    doc.get("Global Run Parameters")->add("load_imbalance", params.load_imbalance);
    if (params.mv_overlap_comm_comp == 1) {
        std::string val("1 (yes)");
        doc.get("Global Run Parameters")->add("mv_overlap_comm_comp", val);
    }
    else {
        std::string val("0 (no)");
        doc.get("Global Run Parameters")->add("mv_overlap_comm_comp", val);
    }
}

void add_configuration_to_yaml(YAML_Doc& doc, int numprocs, int numthreads)
{
    doc.get("Global Run Parameters")->add("number of processors", numprocs);

    doc.add("Platform","");
    doc.get("Platform")->add("hostname","Jack's Laptop");
    doc.get("Platform")->add("kernel name","Irrelevant");
    doc.get("Platform")->add("kernel release","Irrelevant");
    doc.get("Platform")->add("processor","i7");

    doc.add("Build","");
    doc.get("Build")->add("CXX","C++");
#if MINIFE_INFO != 0
    doc.get("Build")->add("compiler version",MINIFE_CXX_VERSION);
#endif
    doc.get("Build")->add("CXXFLAGS","Unimportant");
    std::string using_mpi("no");
#ifdef HAVE_MPI
    using_mpi = "yes";
#endif
    doc.get("Build")->add("using MPI",using_mpi);
}

void add_timestring_to_yaml(YAML_Doc& doc)
{
    std::time_t rawtime;
    struct tm * timeinfo;
    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    std::ostringstream osstr;
    osstr.fill('0');
    osstr << timeinfo->tm_year+1900 << "-";
    osstr.width(2); osstr << timeinfo->tm_mon+1 << "-";
    osstr.width(2); osstr << timeinfo->tm_mday << ", ";
    osstr.width(2); osstr << timeinfo->tm_hour << "-";
    osstr.width(2); osstr << timeinfo->tm_min << "-";
    osstr.width(2); osstr << timeinfo->tm_sec;
    std::string timestring = osstr.str();
    doc.add("Run Date/Time",timestring);
}