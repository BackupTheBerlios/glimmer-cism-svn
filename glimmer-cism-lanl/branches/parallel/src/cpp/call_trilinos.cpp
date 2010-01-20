#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Comm.h"
#include "Epetra_Version.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Trilinos_Util.h"

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"

extern "C" {

  //  void solve_(int& nnz, int& order, int* row, int* col, double* val,
  //	      double* rhs, double* solution, double err, int niters) {
  void solve_(int& nnz, int& order, int* row, int* col, double* val,
	      double* rhs, double* solution) {

    try{
    
#ifdef HAVE_MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    int i, j, ierr, nPEs;
    int MyPID = Comm.MyPID();
    bool verbose = (MyPID == 0);
    Epetra_Map RowMap(order, 0, Comm);
    int NumMyElements = RowMap.NumMyElements();
    int *MyGlobalElements = new int[NumMyElements];
    RowMap.MyGlobalElements(&MyGlobalElements[0]);

    MPI_Comm_size(MPI_COMM_WORLD, &nPEs);

    //-------------------------------------------------------------------------
    // RN_20100120: Counting non-zero entries per row and determining the max
    //-------------------------------------------------------------------------
    int *NumEntriesPerRow = new int[order]; // a cheap solution for now
    int localmax = 0;
    int globalmax = 0;

    for (j=0; j<order; ++j) {
      NumEntriesPerRow[j] = 0;
    }

    // Count non-zero entries per row.
    for (j=0; j<nnz; ++j) {
      if (RowMap.MyGID(row[j]) ) {
	NumEntriesPerRow[row[j] ] += 1;
      }
    }
    
    /*
    // Debug
    cout << "Non-zero counts from " << MyPID << endl;
    for (j=0; j<order; ++j) {
      cout << NumEntriesPerRow[j] << endl;
    }    
    */

    // Determine localmax.
    for (j=0; j<order; ++j) {
      if (NumEntriesPerRow[j] > localmax) {
	localmax = NumEntriesPerRow[j];
      }
    }

    // Determine globalmax.
    if (MyPID == 0) {
      int *MaxFromProcessors = new int[nPEs];
      MPI_Status status;
      MaxFromProcessors[0] = localmax;
      for (j=1; j<nPEs; ++j) {
	// Receive localmax.
	MPI_Recv(&MaxFromProcessors[j], 1, MPI_INT, j, 1, MPI_COMM_WORLD,
		 &status);
      }

      /*
      // Debug
      cout << "Maxes from all the processors:" << endl;
      for (j=0; j<nPEs; ++j) {
	cout << MaxFromProcessors[j] << endl;
      }
      */

      // Find globalmax.
      for (j=0; j<nPEs; ++j) {
	if (MaxFromProcessors[j] > globalmax) {
	  globalmax = MaxFromProcessors[j];
	}
      }
      // Broadcast globalmax.
      MPI_Bcast(&globalmax, 1, MPI_INT, 0, MPI_COMM_WORLD);

      delete[] MaxFromProcessors;
    }
    else {
      // Send localmax.
      MPI_Send(&localmax, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      // Receive globalmax.
      MPI_Bcast(&globalmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    delete[] NumEntriesPerRow;

    //-------------------------------------------------------------------------
    // RN_20091221: Taking care of the matrix
    //-------------------------------------------------------------------------
    //cout << "global max: " << globalmax << endl;

    //Epetra_CrsMatrix A(Copy, RowMap, ColMap, NumEntriesPerRow);
    Epetra_CrsMatrix A(Copy, RowMap, globalmax);

    // Inserting values into the matrix    
    /*
    for (i=0; i<NumMyElements; ++i) {
      for (j=0; j<nnz; ++j) {
	if (MyGlobalElements[i] == row[j]) {
	  ierr = A.InsertGlobalValues(MyGlobalElements[i],
				      1, &(val[j]), &(col[j]));
	  assert(ierr == 0);
	}
      }
    }
    */
    for (j=0; j<nnz; ++j) {
      if (RowMap.MyGID(row[j]) ) {
	ierr = A.InsertGlobalValues(row[j], 1, &(val[j]), &(col[j]) );
	assert(ierr == 0);
      }
    }

    ierr = A.FillComplete();
    assert(ierr == 0);
    
    //-------------------------------------------------------------------------
    // RN_20091221: Taking care of the rhs
    //-------------------------------------------------------------------------
    //Epetra_Vector b(Copy, RowMap, rhs);
    Epetra_Vector b(RowMap);

    // Inserting values into the rhs
    double *MyGlobalValues = new double[NumMyElements];
    for (j=0; j<NumMyElements; ++j) {
      MyGlobalValues[j] = rhs[MyGlobalElements[j] ];
    }
    ierr = b.ReplaceGlobalValues(NumMyElements, &MyGlobalValues[0],
				 &MyGlobalElements[0]);

    //-------------------------------------------------------------------------
    // RN_20091221: Taking care of the solution
    //-------------------------------------------------------------------------
    Epetra_Vector x(RowMap);

    Teuchos::ParameterList paramList;

    Teuchos::RCP<Teuchos::ParameterList>
      paramList1 = Teuchos::rcp(&paramList, false);
    Teuchos::updateParametersFromXmlFile("./strat1.xml", paramList1.get() );

    // For debugging =)
    //cout << "A: " << A << endl;
    //    cout << "b: " << b << endl;
    //    cout << "x: " << x << endl;

    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();

    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;

    Teuchos::RCP<Epetra_CrsMatrix> epetraOper = Teuchos::rcp(&A, false);
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
    //    cout << "A: " << A << endl;
    //    cout << "b: " << b << endl;
    //    cout << "x: " << x << endl;

    //    Epetra_Vector temp(RowMap);
    //    double residualNorm;
    //    A.Multiply(false, x, temp);
    //    temp.Update(-1, b, 1);
    //    temp.Norm2(&residualNorm);

    //    cout << "Residual Norm: " << residualNorm << endl;

    Epetra_LocalMap localMap(order, 0, Comm);
    Epetra_Vector xExtra(localMap); // local vector in each processor

    // RN_20091218: Create an import map and then import the data to the
    // local vector.
    Epetra_Import import(localMap, RowMap);

    xExtra.Import(x, import, Add);
    xExtra.ExtractCopy(solution);

    delete[] MyGlobalElements;
    delete[] MyGlobalValues;
    }
  
    catch (std::exception& e) {
      cout << e.what() << endl;
    }
    catch (string& s) {
      cout << s << endl;
    }
    catch (char *s) {
      cout << s << endl;
    }
    catch (...) {
      cout << "Caught unknown exception!" << endl;
    }
  }

} // extern "C"
