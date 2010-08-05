#include <iostream>
#include "matrixInterface.hpp"

// Constructor
TrilinosMatrix_Interface::TrilinosMatrix_Interface(int bandwidth, int mySize, int* myIndices, const Epetra_Comm& comm)
  : bandwidth_(bandwidth), matrixOrder_(-1), comm_(comm) {
  
  rowMap_ =  Teuchos::rcp(new Epetra_Map(-1, mySize, myIndices, 1, comm) );
  matrixOrder_ = rowMap_->NumGlobalElements();

  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap_, bandwidth) );
  isFillCompleted_ = 0;

  // create map of full vector
  fullMap_ =  Teuchos::rcp(new Epetra_LocalMap(matrixOrder_, 1, comm));
}

// Destructor
TrilinosMatrix_Interface::~TrilinosMatrix_Interface() {
}

// Fix the sparsity patter by calling FillComplete
void TrilinosMatrix_Interface::finalizeSparsity() {
  isFillCompleted_ = true;
  int ierr = operator_->FillComplete();
  assert (ierr==0);
}

// Update the value of bandwidth_ // RN_20100121: probably not needed
void TrilinosMatrix_Interface::updateBandwidth(int bandwidth) {
  bandwidth_ = bandwidth;
}

// Update the operator and also the corresponding row map.
void TrilinosMatrix_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  rowMap_ = Teuchos::rcp(new Epetra_Map(operator_->RowMap() ) );
  isFillCompleted_ = false;
}
