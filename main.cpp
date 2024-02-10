#include "Solver.h"

int main() {

    Solver2D doubleMachReflection;
    doubleMachReflection.setDomainData(4.0, 100, 1.0, 25, 0.2, 0.3);
    doubleMachReflection.doAnalysis("/home/jeeemieee/simulationOutput/doubleMachReflection/testAnalysis/",10);

    return 0;
}
