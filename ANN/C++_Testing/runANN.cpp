#include <iostream>
#include "../onnx_inference/ANN.H"

int main(){

    //trained model information - Mention full path is available
    char file_path[] = "C12_N2_TML_PR.onnx";
    char mins_path[] = "mins_data.txt";
    char maxs_path[] = "maxs_data.txt";

    //Creating the ANN object
    ANN network(file_path,mins_path,maxs_path,4,3);

    //sample input
    std::vector<float> input = {85318.03767507989,4500000.0,0.9801010101010101,0.01989898989898986};

    std::vector<float> output = network.ANN_predict(input);
    std::cout << "Printing \n" << output.size() << "\n";
    

    for(int i=0;i<output.size();i++){
        std::cout << output[i] << "\n";
    }

    return 0;
}