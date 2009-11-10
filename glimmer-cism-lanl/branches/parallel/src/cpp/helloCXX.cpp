#include <iostream>
//#include "Epetra_ConfigDefs.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

//using namespace std;

extern "C"
{
    void hellocxx_ () ;
}

void hellocxx_ () {

	//cout << "  Not doing anything in helloCXX" << endl ;
	cout << "Erin's test" << endl ;
	cout << Epetra_Version() << endl << endl ;

	Epetra_SerialComm Comm;

	int NumElems = 10 ;

	Epetra_Map Map(NumElems, 0, Comm) ;

	Epetra_Vector x(Map) ;
	Epetra_Vector b(Map) ;

	b.Random();
	x.Update(2.0, b, 0.0) ; //x = 2*b

	double bnorm, xnorm ;
	x.Norm2(&xnorm) ;
	b.Norm2(&bnorm) ;

	cout << "2 norm of x = " << xnorm << endl 
	     << "2 norm of b = " << bnorm << endl ;

	return ;

}
