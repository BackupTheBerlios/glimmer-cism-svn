#include <iostream>
#include "Simple_Interface.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

// Define global variables.
static Teuchos::RCP<Simple_Interface> interface;

extern "C" {

  //================================================================
  // RN_20091118: This needs to be called only once in the beginning
  // to set up the problem.
  //================================================================
  void initialize_(int& bandwidth, int& nnz, int& size, int* rowInd, int* colInd, double* val) {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    cout << " ======================================" << endl;
    cout << " IN INITIALIZE()" << endl;

    int i, j, ierr;
    //int bandwidth = 5; // approximately number of entries per row

    // Explicitly change from 1-based indexing to 0-based indexing
    for (i=0; i<nnz; ++i) {
      rowInd[i] = rowInd[i] - 1;
      colInd[i] = colInd[i] - 1;
    }

    interface = Teuchos::rcp(new Simple_Interface(bandwidth, size, comm) );

    const Epetra_Map& map = interface->getMap();
    Teuchos::RCP<Epetra_CrsMatrix> matrix = interface->getOperator();

    int numMyElements = map.NumMyElements();
    int *myGlobalElements = new int[numMyElements];
    map.MyGlobalElements(&myGlobalElements[0]);

    // Inserting values into the matrix    
    for (i=0; i<numMyElements; ++i) {
      for (j=0; j<nnz; ++j) {
	if (myGlobalElements[i] == rowInd[j]) {
	  ierr = matrix->InsertGlobalValues(myGlobalElements[i],
					   1, &(val[j]), &(colInd[j]));
	  assert(ierr == 0);
	}
      }
    }

    ierr = matrix->FillComplete();
    assert(ierr == 0);
    //cout << "Matrix: " << endl << *matrix << endl;

    cout << " ======================================" << endl;
  }

  //============================================================
  // RN_20091118: This is to update the matrix with new entries.
  //============================================================
  void update_(int& rowInd, int& colInd, double& val) {
    cout << " ======================================" << endl;
    cout << " IN UPDATE()" << endl;

    int i, j, ierr = 0;

    // Change from 1-based indexing to 0-based indexing
    rowInd = rowInd - 1;
    colInd = colInd - 1;
    //cout << "A triplet: " << rowInd << " " << colInd << " " << val << endl;

    Teuchos::RCP<Epetra_CrsMatrix> matrix = interface->getOperator();
    ierr = matrix->ReplaceGlobalValues(rowInd, 1, &val, &colInd);
    //    assert(ierr == 0);

    if (ierr != 0) {
      cout << "This new entry (" << rowInd << ", " << colInd << ", "
	   << val << ") did not exist before. A new matrix will be formed!"
	   << endl;

      int matrixSize = interface->matrixOrder();
      const int bandwidth = interface->bandwidth();

#ifdef HAVE_MPI
      Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
      Epetra_SerialComm comm;
#endif

      Epetra_Map map(matrixSize, 0, comm );
      int numMyElements = map.NumMyElements();
      int *myGlobalElements = new int[numMyElements];
      map.MyGlobalElements(&myGlobalElements[0]);

      Teuchos::RCP<Epetra_CrsMatrix> newMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, bandwidth) );

      //const int length = 10;
      int numEntries;
      double *values = new double[bandwidth];
      int *indices = new int[bandwidth];

      // Copy the old matrix to the new matrix.
      for (i=0; i<numMyElements; ++i) {
	for (j=0; j<matrixSize; ++j) {
	  if (myGlobalElements[i] == j) {
	    ierr = matrix->ExtractGlobalRowCopy(j, bandwidth, numEntries, values, indices);
	    assert(ierr == 0);
	    ierr = newMatrix->InsertGlobalValues(myGlobalElements[i],
						numEntries, &(values[0]),
						&(indices[0]) );
	    assert(ierr == 0);
	  }
	}

      }

      // Insert the new entry.
      for (i=0; i<numMyElements; ++i) {
	if (myGlobalElements[i] == rowInd) {
	  ierr = newMatrix->InsertGlobalValues(rowInd, 1, &val, &colInd);
	  assert(ierr == 0);
	}
      }

      ierr = newMatrix->FillComplete();
      assert(ierr == 0);
      //cout << "New matrix: " << endl << newMatrix << endl;

      interface->updateOperator(newMatrix);
      //Teuchos::RCP<Epetra_CrsMatrix> temp = interface->getOperator();
      //cout << "Updated matrix: " << endl << *temp << endl;

      delete[] values;
      delete[] indices;

    cout << " ======================================" << endl;
    }
  }

  //========================================================
  // RN_20091118: This is to make calls to Trilinos solvers.
  //========================================================
  void solve_(double* rhs, double* solution) {
    cout << " ======================================" << endl;
    cout << " IN SOLVE()" << endl;

    Teuchos::RCP<Epetra_CrsMatrix> epetraOper = interface->getOperator();
    const Epetra_Map& map = epetraOper->RowMap();
    //const Epetra_Map& map = interface->getMap();

    Epetra_Vector b(Copy, map, rhs);
    Epetra_Vector x(map);

    //cout << "rhs" << endl << b << endl;

    Teuchos::ParameterList paramList;

    Teuchos::RCP<Teuchos::ParameterList>
      paramList1 = Teuchos::rcp(&paramList, false);
    Teuchos::updateParametersFromXmlFile("strat1.xml", paramList1.get() );

    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    Teuchos::RCP<Epetra_Vector> epetraRhs = Teuchos::rcp(&b, false);
    Teuchos::RCP<Epetra_Vector> epetraSol = Teuchos::rcp(&x, false);

    Teuchos::RCP<const Thyra::LinearOpBase<double> >
      thyraOper = Thyra::epetraLinearOp(epetraOper);
    Teuchos::RCP<Thyra::VectorBase<double> >
      thyraRhs = Thyra::create_Vector(epetraRhs, thyraOper->range() );
    Teuchos::RCP<Thyra::VectorBase<double> >
      thyraSol = Thyra::create_Vector(epetraSol, thyraOper->domain() );

    linearSolverBuilder.setParameterList(Teuchos::rcp(&paramList, false) );

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
      lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");

    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> >
      lows = Thyra::linearOpWithSolve(*lowsFactory, thyraOper);

    Thyra::SolveStatus<double>
      status = Thyra::solve(*lows, Thyra::NOTRANS, *thyraRhs, &*thyraSol);

    thyraSol = Teuchos::null;

    // For debugging =)
    //    cout << "A: " << *epetraOper << endl;
    //    cout << "b: " << b << endl;
    //    cout << "x: " << x << endl;

    //    Epetra_Vector temp(rangeMap);
    //    double residualNorm;
    //    matrix.Multiply(false, x, temp);
    //    temp.Update(-1, b, 1);
    //    temp.Norm2(&residualNorm);

    //    cout << "Residual Norm: " << residualNorm << endl;

    x.ExtractCopy(solution);

    cout << " ======================================" << endl;
  }

  //======================
  // RN_20091118: Wrap up!
  //======================
  void finalize_() {
    cout << " ======================================" << endl;
    cout << " IN FINALIZE()" << endl;

    // Might want to shift back to 1-based indexing?

    cout << " ======================================" << endl;
  }

} // extern"C"
