#include <iostream>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
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

class Simple_Interface {
public:
  // Constructor
  Simple_Interface(int bandwidth, int matrixOrder, const Epetra_Comm& comm);

  // Destructor
  ~Simple_Interface();

  // Accessors
  int isFillCompleted() {return isFillCompleted_;};
  const int bandwidth() const {return bandwidth_;};
  const int matrixOrder() const {return matrixOrder_;};
  const Epetra_Map& getMap() const {return *globalMap_;};
  Teuchos::RCP<Epetra_CrsMatrix>& getOperator() {return operator_;};

  // Mutators
  void updateFill(int fill);
  void updateBandwidth(int bandwidth); // RN_20100121: probably not needed
  void updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator);

private:
  int isFillCompleted_; // to indicate if operator_ is "FillComplete()"ed
  int bandwidth_;
  int matrixOrder_;
  const Epetra_Comm& comm_;
  Teuchos::RCP<Epetra_CrsMatrix> operator_;
  Teuchos::RCP<Epetra_Map> globalMap_;
};
