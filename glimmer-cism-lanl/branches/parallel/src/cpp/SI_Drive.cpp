#include <iostream>
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_CombineMode.h"
#include "Simple_Interface.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_Time.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

// Define global variables.
static Teuchos::RCP<Simple_Interface> interface;

extern "C" {

  //================================================================
  // RN_20091215: This needs to be called only once in the beginning
  // to set up the problem.
  //================================================================
  void initialize_(int& bandwidth, int& size) {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    cout << " ======================================" << endl;
    cout << " IN INITIALIZE()" << endl;

    // Create an interface that holds a CrsMatrix instance.
    interface = Teuchos::rcp(new Simple_Interface(bandwidth, size, comm) );

    cout << " ======================================" << endl;
  }

  //============================================================
  // RN_20100201: This is to check if this entry already exists.
  //============================================================
  void exist_(int& rowInd, int& colInd, int& flag) {
    // RN_20100201: This is not needed since the sparsity pattern
    // does not change within a Picard iteration.
  }

  //============================================================
  // RN_20091118: This is to update the matrix with new entries.
  //============================================================
  void update_(int& rowInd, int& colInd, double& val) {
    //cout << " ======================================" << endl;
    //cout << " IN UPDATE()" << rowInd << endl;

    /*    
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    */

    //int i;
    int j, ierr = 0;
    int rowInd1 = rowInd - 1;
    int colInd1 = colInd - 1;

    // Change from 1-based indexing to 0-based indexing
    //cout << "A triplet 1: " << rowInd << " " << colInd << " " << val << endl;
    //cout << "A triplet 2: " << rowInd1 << " " << colInd1 << " " << val << endl;

    int fill = interface->isFillCompleted();
    Teuchos::RCP<Epetra_CrsMatrix> matrix = interface->getOperator();
    Epetra_Map map = interface->getMap();
    //cout << "GID: " << map.GID(2) << ", and LID: " << map.LID(2) << endl;
    //cout << "And " << map.MyGID(2) << endl;
    int numMyElements = map.NumMyElements();
    //int *myGlobalElements = new int[numMyElements];
    //map.MyGlobalElements(&myGlobalElements[0]);

    /*
    for (i=0; i<numMyElements; ++i) {
      if (myGlobalElements[i] == rowInd1) {
	ierr = matrix->ReplaceGlobalValues(rowInd1, 1, &val, &colInd1);
	//assert(ierr == 0);
	break;
      }
    }
    */
    int pidList, lidList, junk;
    //junk = map.RemoteIDList(1, &rowInd1, &pidList, &lidList);
    if (map.MyGID(rowInd1) ) {
      ierr = matrix->ReplaceGlobalValues(rowInd1, 1, &val, &colInd1);
    }
    // Is there any way to improve this?
    //    cout << "MPI_Broadcast: " << ierr << ", "<< pidList << ", "
    //	 << rowInd1 << endl;
    
    //MPI_Bcast(&ierr, 1, MPI_INT, pidList, MPI_COMM_WORLD);
    //cout << "ierr" << ierr << endl;
    //cout << "why did you stop?" << endl;

    
    if (ierr != 0) { // Sparsity pattern has changed.
      //cout << "i don't know" << endl;
      if (fill == 0) { // The matrix has not been "FillComplete()"ed.
	//	cout << "This new entry (" << rowInd1 << ", " << colInd1 << ", "
	//   << val << ") did not exist before. No new matrix is formed!"
	//   << endl;

	if (map.MyGID(rowInd1) ) {
	  ierr = matrix->InsertGlobalValues(rowInd1, 1, &val, &colInd1);
	}
      }
      else { // The matrix is "FillComplete()"ed. A new matrix is needed.      
	cout << "This new entry (" << rowInd1 << ", " << colInd1 << ", "
	     << val << ") did not exist before. A new matrix will be formed!"
	     << endl;

	int matrixSize = interface->matrixOrder();
	int bandwidth = interface->bandwidth();
	
	Teuchos::RCP<Epetra_CrsMatrix> newMatrix =
	  Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, bandwidth) );
	
	//const int length = 10;
	int numEntries;
	double *values = new double[bandwidth];
	int *indices = new int[bandwidth];
	
	// Copy the old matrix to the new matrix.
	for (j=0; j<matrixSize; ++j) {
	  if (map.MyGID(j) ) {
	    int aNumber = 20;
	    ierr = matrix->ExtractGlobalRowCopy(j, aNumber, numEntries,
						values, indices);
	    //ierr = matrix->ExtractGlobalRowCopy(j, bandwidth, numEntries,
	    //				values, indices);
	    cout << ierr << endl;
	    assert(ierr >= 0);
	    ierr = newMatrix->InsertGlobalValues(j, numEntries, &(values[0]),
						 &(indices[0]) );
	    cout << ierr << endl;
	    assert(ierr >= 0);
	  }
	}
	
	// Insert the new entry.
	if (map.MyGID(rowInd1) ) {
	  ierr = newMatrix->InsertGlobalValues(rowInd1, 1, &val, &colInd1);
	}

	//ierr = newMatrix->FillComplete();
	//assert(ierr == 0);
	//cout << "New matrix: " << endl << newMatrix << endl;
	
	interface->updateOperator(newMatrix);
	interface->updateFill(0);
	
	//Teuchos::RCP<Epetra_CrsMatrix> temp = interface->getOperator();
	//cout << "Updated matrix: " << endl << *temp << endl;
	
	delete[] values;
	delete[] indices;
      } // else
    
    }
    
    //cout << " ======================================end" << endl;
  }

  //========================================================
  // RN_20091118: This is to make calls to Trilinos solvers.
  //========================================================
  void differentsolve_(double* rhs, double* solution, double& elapsedTime) {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    cout << " ======================================" << endl;
    cout << " IN SOLVE()" << endl;

    // RN_20100211: Start timing
    Teuchos::Time linearTime("LinearTime");
    linearTime.start();

    int j, ierr;
    Teuchos::RCP<Epetra_CrsMatrix> epetraOper = interface->getOperator();

    // For debugging =)
    //cout << "A: " << *epetraOper << endl;

    int fill = interface->isFillCompleted();
    if (fill == 0) {
      //cout << "fill" << fill << endl;
      ierr = epetraOper->FillComplete();
      assert(ierr == 0);
      interface->updateFill(1);
      //fill = interface->isFillCompleted();
      //cout << "fill" << fill << endl;
    }

    // For debugging =)
    //cout << "A: " << *epetraOper << endl;

    const Epetra_Map& map = epetraOper->RowMap();
    //const Epetra_Map& map = interface->getMap();
    int numMyElements = map.NumMyElements();
    int *myGlobalElements = new int[numMyElements];
    map.MyGlobalElements(&myGlobalElements[0]);

    //Epetra_Vector b(Copy, map, rhs);
    Epetra_Vector b(map);

    // Inserting values into the rhs.
    double *myGlobalValues = new double[numMyElements];
    for (j=0; j<numMyElements; ++j) {
      myGlobalValues[j] = rhs[myGlobalElements[j] ];
    }
    ierr = b.ReplaceGlobalValues(numMyElements, &myGlobalValues[0],
				 &myGlobalElements[0]);

    // For debugging =)
    //cout << "b: " << b << endl;

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
    //cout << "x: " << x << endl;

    //    Epetra_Vector temp(rangeMap);
    //    double residualNorm;
    //    matrix.Multiply(false, x, temp);
    //    temp.Update(-1, b, 1);
    //    temp.Norm2(&residualNorm);

    //    cout << "Residual Norm: " << residualNorm << endl;

    int matrixSize = interface->matrixOrder();
    Epetra_LocalMap localMap(matrixSize, 0, comm);
    Epetra_Vector xExtra(localMap); // local vector in each processor

    // RN_20091218: Create an import map and then import the data to the
    // local vector.
    Epetra_Import import(localMap, map);

    xExtra.Import(x, import, Add);
    xExtra.ExtractCopy(solution);

    delete[] myGlobalElements;

    elapsedTime = linearTime.stop();
    cout << "Total time elapsed for calling Solve(): " << elapsedTime << endl;

    cout << " ======================================" << endl;
  }

  //======================
  // RN_20091118: Wrap up!
  //======================
  void finalize_() {
    cout << " ======================================" << endl;
    cout << " IN FINALIZE()" << endl;

    // The following is not essential, but only for completion.
    int ierr;
    Teuchos::RCP<Epetra_CrsMatrix> epetraOper = interface->getOperator();
    int fill = interface->isFillCompleted();
    if (fill == 0) {
      //cout << "I am in here!" << endl;
      ierr = epetraOper->FillComplete();
      assert(ierr == 0);
      interface->updateFill(1);
      fill = interface->isFillCompleted();
    }

    //cout << "Matrix: " << endl << *epetraOper << endl;

    cout << " ======================================" << endl;
  }

} // extern"C"
