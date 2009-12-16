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
    
#ifdef EPETRA_MPI
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    
    int i, j, ierr, maxID, max;
    //    int MyPID = Comm.MyPID();
    //    bool verbose = (MyPID == 0);
    Epetra_Map RowMap(order, 0, Comm);
    int NumMyElements = RowMap.NumMyElements();
    int *MyGlobalElements = new int[NumMyElements];
    RowMap.MyGlobalElements(&MyGlobalElements[0]);
    Epetra_Map ColMap(order, 0, Comm);
    ColMap.MyGlobalElements(&MyGlobalElements[0]);

    // E>> trying to make NumEntriesPerRow array
    int *NumEntriesPerRow = new int[NumMyElements];
    for (i=0; i<NumMyElements; i++) { NumEntriesPerRow[i] = 0;}
    //NumEntriesPerRow = 0;
    maxID= 0; max = 1;
    for (i=0; i<nnz; i++) {	
	if ( row[i] >= 0) {
	    NumEntriesPerRow[row[i]] += 1;
	    if ( NumEntriesPerRow[row[i]] > max) max = NumEntriesPerRow[row[i]];
	    if ( row[i] > maxID) maxID = row[i];
	}
    }

    maxID = 0;
    for (i=0; i<nnz; i++) {	
	    if ( col[i] > maxID) maxID = col[i];
    }
    delete[] NumEntriesPerRow;

    // the number of entries per row in the matrix
    int two = 2;
    //Epetra_CrsMatrix A(Copy, RowMap, ColMap, NumEntriesPerRow);
    Epetra_CrsMatrix A(Copy, RowMap, max);
    Epetra_Vector b(Copy, RowMap, rhs);
    Epetra_Vector x(RowMap);

    Teuchos::ParameterList paramList;

    Teuchos::RCP<Teuchos::ParameterList>
      paramList1 = Teuchos::rcp(&paramList, false);
    Teuchos::updateParametersFromXmlFile("strat1.xml", paramList1.get() );

    // Inserting values into the matrix    
    for (i=0; i<NumMyElements; ++i) {
      for (j=0; j<nnz; ++j) {
	if (MyGlobalElements[i] == row[j]) {
	  ierr = A.InsertGlobalValues(MyGlobalElements[i],
				      1, &(val[j]), &(col[j]));
	  assert(ierr == 0);
	}
      }
    }

    ierr = A.FillComplete();
    assert(ierr == 0);
    
    // For debugging =)
    //    cout << "A: " << A << endl;
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

    x.ExtractCopy(solution);
  }

} // extern "C"
