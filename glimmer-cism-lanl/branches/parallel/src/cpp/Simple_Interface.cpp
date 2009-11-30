#include <iostream>
#include "Simple_Interface.hpp"

// Constructor
Simple_Interface::Simple_Interface(int bandwidth, int matrixSize, const Epetra_Comm& comm) : bandwidth_(bandwidth), matrixOrder_(matrixSize), comm_(comm) {
  globalMap_ =  Teuchos::rcp(new Epetra_Map(matrixSize, 0, comm) );
  
  operator_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalMap_, bandwidth) );
}

// Destructor
Simple_Interface::~Simple_Interface() {
}

// Update the operator and also the corresponding global map.
void Simple_Interface::updateOperator(Teuchos::RCP<Epetra_CrsMatrix> newOperator) {
  operator_ = newOperator;
  globalMap_ = Teuchos::rcp(new Epetra_Map(operator_->RowMap() ) );
}
