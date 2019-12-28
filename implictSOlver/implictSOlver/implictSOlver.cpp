#include "stdafx.h"
#include "temSolver.h"

implictSolver solver;



int _tmain(int argc, _TCHAR* argv[])
{

	solver.set_init_value();

	solver.allo_mem();

	solver.init_field();

	solver.run();

	return 0;
}

