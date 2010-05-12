#include <iostream>
#include "Simple_Interface.hpp"

// Constructor (NO LONGER USED)
Simple_Interface::Simple_Interface(int bandwidth, int matrixSize, const Epetra_Comm& comm)
  : bandwidth_(bandwidth), matrixOrder_(matrixSize), comm_(comm) {
  rowMap_ =  Teuchos::rcp(new Epetra_Map(matrixSize, 1, comm) );
  
  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap_, bandwidth) );
  isFillCompleted_ = 0;

  // create map of full vector
  fullMap_ =  Teuchos::rcp(new Epetra_LocalMap(matrixSize, 1, comm));
}

// Constructor
Simple_Interface::Simple_Interface(int bandwidth, int mySize, int* myIndices, const Epetra_Comm& comm)
  : bandwidth_(bandwidth), matrixOrder_(-1), comm_(comm) {
  
  rowMap_ =  Teuchos::rcp(new Epetra_Map(-1, mySize, myIndices, 1, comm) );
  matrixOrder_ = rowMap_->NumGlobalElements();

  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap_, bandwidth) );
  isFillCompleted_ = 0;

  // create map of full vector
  fullMap_ =  Teuchos::rcp(new Epetra_LocalMap(matrixOrder_, 1, comm));
}

// Destructor
Simple_Interface::~Simple_Interface() {
}

// Update the value of fill_
void Simple_Interface::updateFill(int fill) {
  isFillCompleted_ = fill;
}

// Update the value of bandwidth_ // RN_20100121: probably not needed
void Simple_Interface::updateBandwidth(int bandwidth) {
  bandwidth_ = bandwidth;
}

// Update the operator and also the corresponding row map.
void Simple_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  rowMap_ = Teuchos::rcp(new Epetra_Map(operator_->RowMap() ) );
}
