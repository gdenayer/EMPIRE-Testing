In cavityFSIBenchmarkOF/system/fvSolution, it is very important to set nNonOrthogonalCorrectors as 

PIMPLE
{
    nNonOrthogonalCorrectors 2;
}

otherwise, the OF computation terminates at some time with nan error.