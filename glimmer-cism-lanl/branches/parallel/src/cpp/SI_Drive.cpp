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
static Teuchos::RCP<Epetra_Map> partitionMap;

extern "C" {

  // doPartition and getPartition use Epetra to partition the global
  // problem in parallel. These functions are temporary as the partition
  // will be chosen in glimmer, but allow us to mature the initTrilinos
  // interface.
  void dopartition_(int& matrixSize, int& mySize) {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    partitionMap = Teuchos::rcp(new Epetra_Map(matrixSize,1,comm) );
    mySize = partitionMap->NumMyElements();

    cout << "XXX doPartition has mySize = " << mySize << endl;
  }

  void getpartition_(int& mySize, int* myIndicies) {

      // Copy indices into array to send back to glimmer
      partitionMap->MyGlobalElements(myIndicies);

  }
  //================================================================
  //================================================================
  // RN_20091215: This needs to be called only once in the beginning
  // to set up the problem.
  //================================================================
  void inittrilinos_(int& bandwidth, int& mySize, int* myIndicies) {
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    cout << " ======================================" << endl;
    cout << " IN INITIALIZE()" << endl;

    // Create an interface that holds a CrsMatrix instance.
    interface = Teuchos::rcp(new Simple_Interface(bandwidth, mySize, myIndicies, comm) );

    cout << " ======================================" << endl;
  }

  // initialize can be removed soon -- no longer used
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

    int fill = interface->isFillCompleted();
    Teuchos::RCP<Epetra_CrsMatrix>& matrix = interface->getOperator();
    const Epetra_Map& map = matrix->RowMap();
    int numMyElements = map.NumMyElements();

    int pidList, lidList, junk;
    if (map.MyGID(rowInd) ) {
      ierr = matrix->ReplaceGlobalValues(rowInd, 1, &val, &colInd);
    }
    
    
    if (ierr != 0) { // Sparsity pattern has changed.
      if (fill == 0) { // The matrix has not been "FillComplete()"ed.
	if (map.MyGID(rowInd) ) {
	  ierr = matrix->InsertGlobalValues(rowInd, 1, &val, &colInd);
	}
      }
      else { // The matrix is "FillComplete()"ed. A new matrix is needed.      
	cout << "This new entry (" << rowInd << ", " << colInd << ", "
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
	    int aNumber = bandwidth;
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
	if (map.MyGID(rowInd) ) {
	  ierr = newMatrix->InsertGlobalValues(rowInd, 1, &val, &colInd);
	}

	interface->updateOperator(newMatrix);
	interface->updateFill(0);
	
	delete[] values;
	delete[] indices;
      }
    
    }
    
  }

  //========================================================
  // RN_20091118: This is to make calls to Trilinos solvers.
  //========================================================
  void differentsolve_(double* rhs, double* answer, double& elapsedTime) {
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
    Teuchos::RCP<Epetra_CrsMatrix>& epetraOper = interface->getOperator();

    int fill = interface->isFillCompleted();
    if (fill == 0) {
      ierr = epetraOper->FillComplete();
      assert(ierr == 0);
      interface->updateFill(1);
    }

    const Epetra_Map& map = epetraOper->RowMap();
    int numMyElements = map.NumMyElements();
    int *myGlobalElements = new int[numMyElements];
    map.MyGlobalElements(&myGlobalElements[0]);

    Epetra_Vector b(map);

    // Inserting values into the rhs: "-1" for Fortran to C Numbering
    double *myGlobalValues = new double[numMyElements];
    for (j=0; j<numMyElements; ++j) {
      myGlobalValues[j] = rhs[myGlobalElements[j] -1 ];
    }
    ierr = b.ReplaceGlobalValues(numMyElements, &myGlobalValues[0],
				 &myGlobalElements[0]);

    Epetra_Vector x(map);

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

#ifdef ONLY_RETURN_OWNED_VECTOR
    x.ExtractCopy(answer);
#else
    // Import result to map with a Halo -- now the global problem
    const Epetra_Map& answerMap = interface->getFullMap();
    Epetra_Vector xExtra(answerMap); // local vector in each processor

    // RN_20091218: Create an import map and then import the data to the
    // local vector.
    Epetra_Import import(answerMap, map);

    xExtra.Import(x, import, Add);
    xExtra.ExtractCopy(answer);
#endif

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
    Teuchos::RCP<Epetra_CrsMatrix>& epetraOper = interface->getOperator();
    int fill = interface->isFillCompleted();
    if (fill == 0) {
      //cout << "I am in here!" << endl;
      ierr = epetraOper->FillComplete();
      assert(ierr == 0);
      interface->updateFill(1);
      fill = interface->isFillCompleted();
    }

    cout << " ======================================" << endl;
  }

} // extern"C"
