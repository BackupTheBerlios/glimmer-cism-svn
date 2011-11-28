#include <iostream>
#include "Epetra_LocalMap.h"
#ifdef HAVE_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif
#include "Teuchos_RCP.hpp"


// Define global variables.
static Teuchos::RCP<const Epetra_Map> partitionMap;

// Jeff's Quick Hash Table GlobalIDs solution.  hash_map is defined in PGI at "include/CC/hash_map"
#include <hash_map>
//map #include <map>

hash_map <int, int> IDs;
// map <int, int> IDs;

extern "C" {

  void globalidsadd_(const int &key, const int &indx) {
    IDs[key] = indx;
  }

  void globalidsclear_(void) {
    IDs.clear();
  }

  int globalidsget_(const int &key) {
    return IDs[key];
  }
// } // extern "C"

// extern "C" {

  // doPartition and getPartition use Epetra to partition the global
  // problem in parallel. 
  // This is needed for serial-glimmer and parallel-trilinos,
  // These function will not be needed with a distributed glimmer
  // since the PDE code will pick the partitioning.

  void dopartition_(int& matrixSize, int& mySize) {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    partitionMap = Teuchos::rcp(new Epetra_Map(matrixSize,1,comm) );
    mySize = partitionMap->NumMyElements();

    cout << "Trilinos Partition: doPartition has mySize = " << mySize << endl;
  }

  void getpartition_(int& mySize, int* myIndicies) {

      // Copy indices into array to send back to glimmer
      partitionMap->MyGlobalElements(myIndicies);

  }
} // extern"C"
