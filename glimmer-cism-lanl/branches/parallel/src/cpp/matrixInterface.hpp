#include <iostream>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"

class TrilinosMatrix_Interface {
public:
  // Constructor
  TrilinosMatrix_Interface(int bandwidth, int mySize, int* myIndices, const Epetra_Comm& comm);

  // Destructor
  ~TrilinosMatrix_Interface();

  // Accessors
  bool isSparsitySet() {return isFillCompleted_;};
  const int bandwidth() const {return bandwidth_;};
  const int matrixOrder() const {return matrixOrder_;};
  const Epetra_Map& getFullMap() const {return *fullMap_;};
  Teuchos::RCP<Epetra_CrsMatrix>& getOperator() {return operator_;};

  // Mutators
  void finalizeSparsity(); // Call FillComplet to lock in sparsity pattern
  void updateBandwidth(int bandwidth); // RN_20100121: probably not needed
  void updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator);

private:
  bool isFillCompleted_; // to indicate if operator_ is "FillComplete()"ed
  int bandwidth_;
  int matrixOrder_;
  const Epetra_Comm& comm_;
  Teuchos::RCP<Epetra_CrsMatrix> operator_;
  Teuchos::RCP<Epetra_Map> rowMap_;
  Teuchos::RCP<Epetra_Map> fullMap_;
};
