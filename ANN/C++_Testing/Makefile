# Define the compiler
CXX = g++

# Define the paths
INCLUDE_PATH = ../onnxruntime-linux-x64-1.14.1/include/
LIB_PATH = ../onnx/onnxruntime-linux-x64-1.14.1/lib/
LIB_NAME = libonnxruntime.so.1.14.1

# Define the source and output files
SRC = runANN.cpp
OUT = runANN

# Compilation rule
$(OUT): $(SRC)
	$(CXX) $(SRC) -I$(INCLUDE_PATH) $(LIB_PATH)$(LIB_NAME) -o $(OUT)

# Clean rule
clean:
	rm -f $(OUT)
