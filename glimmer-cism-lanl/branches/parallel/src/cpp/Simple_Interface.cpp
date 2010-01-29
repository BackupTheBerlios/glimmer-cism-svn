#include <iostream>
#include "Simple_Interface.hpp"

// Constructor
Simple_Interface::Simple_Interface(int bandwidth, int matrixSize, const Epetra_Comm& comm) : bandwidth_(bandwidth), matrixOrder_(matrixSize), comm_(comm) {
  globalMap_ =  Teuchos::rcp(new Epetra_Map(matrixSize, 0, comm) );
  
  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalMap_, bandwidth) );
  isFillCompleted_ = 0;
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

// Update the operator and also the corresponding global map.
void Simple_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  globalMap_ = Teuchos::rcp(new Epetra_Map(operator_->RowMap() ) );
}
