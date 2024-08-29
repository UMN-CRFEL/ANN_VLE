#include<vector>

//constructors
ANN::ANN(char* file_path){
    std::string file(file_path);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    session = Ort::Session(env, file_path, session_options);
}

ANN::ANN(){
}

ANN::ANN(char* file_path, int num_inputs, int num_outputs){
    std::string file(file_path);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    session = Ort::Session(env, file_path, session_options);
    input_num = (size_t)num_inputs;
    output_num = (size_t)num_outputs;
}

ANN::ANN(char* file_path, char* mins_path, char* maxs_path, int num_inputs, int num_outputs){
    std::string file(file_path);
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
    session = Ort::Session(env, file_path, session_options);
    input_num = (size_t)num_inputs;
    output_num = (size_t)num_outputs;

    readFiles(mins_path, 0);
    readFiles(maxs_path, 1);
}

// void ANN::operator=(ANN& s){
//     // session(s.session);
//     // input_num = s.input_num;
//     // output_num = s.output_num;
// }

void ANN::readFiles(char* file_path, int flag){
    //flag ==0 -> mins, flag == 1 -> maxs

    std::string filename(file_path);
    std::ifstream inputFile(filename);

    //counting the number of entries in the file
    int lineCount = 0;
    std::string line;
    while (std::getline(inputFile, line)) {
        ++lineCount;
    }

    int total_num = lineCount;

    if(flag == 0)
    {
        mins = new float[total_num];
        std::ifstream inputFile(filename);
        std::cout<<"Reading Mins from "<<filename<<"\n";

        if (inputFile.is_open()) 
        {
            float value;
            int index = 0;
            
            // Read floats from the file and store them in the array
            while (inputFile >> value && index < total_num) 
            {
                mins[index] = value;
                index++;
            }
        }
        
        // Close the file
        inputFile.close();

        for(int i = 0; i<total_num;i++){
            std::cout<<mins[i]<<"\t";
        }
        std::cout<<std::endl;
    
    }
    else if (flag == 1)
    {
        maxs = new float[total_num];

        std::cout<<"Reading Maxs \n";
        std::ifstream inputFile(filename);

        if (inputFile.is_open()) 
        {
            float value;
            int index = 0;
            
            // Read floats from the file and store them in the array
            while (inputFile >> value && index < total_num) 
            {
                maxs[index] = value;
                index++;
            }
        }
        
        // Close the file
        inputFile.close();

        for (int i = 0; i<total_num;i++){
            std::cout<<maxs[i]<<"\t";
        }
        std::cout<<std::endl;
    }
}

//ANN predict function
std::vector<float> ANN::ANN_predict(std::vector<float>& inputValues){
    //Memory Info
    auto memoryinfo = Ort::MemoryInfo::CreateCpu(OrtDeviceAllocator,OrtMemTypeCPU); 
    const char* inputNames[] = {"input"};
    const char* outputNames[] = {"outputs"};

    auto inputshape = session.GetInputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
    int batchsize = 1;
    inputshape[0] = batchsize;

    std::vector<float>inputs = preprocess(inputValues);
    auto inputOnnxTensor = Ort::Value::CreateTensor<float>(memoryinfo,inputs.data(),inputs.size(),inputshape.data(),inputshape.size());
    auto outputValues = session.Run(Ort::RunOptions(nullptr),inputNames,&inputOnnxTensor,1,outputNames,1);
    // std::cout<<"Inference Complete \n";
    
    float* float_output = outputValues.front().GetTensorMutableData<float>();

    std::vector<float> outputVals;

    for(int i=0; i<output_num; i++){
        outputVals.push_back(float_output[i]);
    }

    std::vector<float> outputs = postprocess(outputVals);
    
    return outputs;
}

std::vector<float> ANN::preprocess(std::vector<float>& inputValues){
    std::vector<float> Pp_inputs;
    float val = 0;
    for(int i = 0; i < input_num; ++i ){
        val = (inputValues[i] - mins[i] )/ (maxs[i] - mins[i]);
        Pp_inputs.push_back(val);
    }
    return Pp_inputs;
}

std::vector<float> ANN::postprocess(std::vector<float>& outputValues){
    std::vector<float> Pp_outputs;
    float val = 0;
    for(int i = 0; i < output_num; ++i ){
        val = outputValues[i]  * (maxs[i+output_num+1] - mins[i+output_num+1]) + mins[i+output_num+1];
        Pp_outputs.push_back(val);
    }
    return Pp_outputs;
}