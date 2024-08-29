#include <iostream>
#include <fstream>
#include "../onnx/onnxruntime-linux-x64-1.14.1/include/onnxruntime_cxx_api.h"
#include <vector>

class ANN{

    private:
    const char* path;
    static const auto inputshape = 0, outputshape = 0;
    size_t input_num, output_num;

    float* mins; 
    float* maxs;

    Ort::Env env;
    Ort::Session session{nullptr};
    Ort::SessionOptions session_options;
    Ort::AllocatorWithDefaultOptions allocator;
    // static const auto memoryinfo;


    Ort::Value input_tensor_{nullptr};
    Ort::Value output_tensor_{nullptr};

    public:
    inline ANN();
    inline ANN(char* file_path);
    inline ANN(char* file_path, int num_inputs, int num_outputs);
    inline ANN(char* file_path, char* mins_path, char* maxs_path, int num_inputs, int num_outputs);
    inline std::vector<float> ANN_predict(std::vector<float>& inputValues);
    inline std::vector<float> preprocess(std::vector<float>& inputValues);
    inline std::vector<float> postprocess(std::vector<float>& outputValues);
    inline void readFiles(char* file_path, int flag);
    // void operator=(ANN& s);
};

#include "ANN.C"

