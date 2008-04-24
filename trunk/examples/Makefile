include ../makevars

EXS = AnalyticFunction LocalOptimization1D NumericalDerivativeND CutOptimization FunctionInLine LocalOptimizationND NumericalDerivative1D

all: $(EXS)

%.C:

$(EXS): %: %.o
	$(CXX) $@.o -o $@ 

clean:        
	rm *.o $(EXS) 
