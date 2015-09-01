#include "NNtest.h"
#include <cmath>

double NNtest::Value(int index,double in0,double in1,double in2) {
   input0 = (in0 - 0.694305)/0.0453022;
   input1 = (in1 - 0.0142609)/0.0785499;
   input2 = (in2 - 2.22695)/2.61238;
   switch(index) {
     case 0:
         return neuron0x25ebf90();
     default:
         return 0.;
   }
}

double NNtest::Value(int index, double* input) {
   input0 = (input[0] - 0.694305)/0.0453022;
   input1 = (input[1] - 0.0142609)/0.0785499;
   input2 = (input[2] - 2.22695)/2.61238;
   switch(index) {
     case 0:
         return neuron0x25ebf90();
     default:
         return 0.;
   }
}

double NNtest::neuron0x25d9db0() {
   return input0;
}

double NNtest::neuron0x25da0f0() {
   return input1;
}

double NNtest::neuron0x25eb080() {
   return input2;
}

double NNtest::input0x25eb460() {
   double input = -502.692;
   input += synapse0x25da4c0();
   input += synapse0x25eb710();
   input += synapse0x25eb750();
   return input;
}

double NNtest::neuron0x25eb460() {
   double input = input0x25eb460();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double NNtest::input0x25eb790() {
   double input = -3.33545;
   input += synapse0x25ebad0();
   input += synapse0x25ebb10();
   input += synapse0x25ebb50();
   return input;
}

double NNtest::neuron0x25eb790() {
   double input = input0x25eb790();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double NNtest::input0x25ebb90() {
   double input = 417.014;
   input += synapse0x25ebed0();
   input += synapse0x25ebf10();
   input += synapse0x25ebf50();
   return input;
}

double NNtest::neuron0x25ebb90() {
   double input = input0x25ebb90();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double NNtest::input0x25ebf90() {
   double input = -42.7354;
   input += synapse0x25ec1b0();
   input += synapse0x25ec1f0();
   input += synapse0x25ec230();
   return input;
}

double NNtest::neuron0x25ebf90() {
   double input = input0x25ebf90();
   return (input * 1)+0;
}

double NNtest::synapse0x25da4c0() {
   return (neuron0x25d9db0()*164.774);
}

double NNtest::synapse0x25eb710() {
   return (neuron0x25da0f0()*-407.048);
}

double NNtest::synapse0x25eb750() {
   return (neuron0x25eb080()*-728.504);
}

double NNtest::synapse0x25ebad0() {
   return (neuron0x25d9db0()*0.106345);
}

double NNtest::synapse0x25ebb10() {
   return (neuron0x25da0f0()*-1.04469);
}

double NNtest::synapse0x25ebb50() {
   return (neuron0x25eb080()*-2.18794);
}

double NNtest::synapse0x25ebed0() {
   return (neuron0x25d9db0()*-0.066145);
}

double NNtest::synapse0x25ebf10() {
   return (neuron0x25da0f0()*166.709);
}

double NNtest::synapse0x25ebf50() {
   return (neuron0x25eb080()*448.754);
}

double NNtest::synapse0x25ec1b0() {
   return (neuron0x25eb460()*-0.0333819);
}

double NNtest::synapse0x25ec1f0() {
   return (neuron0x25eb790()*3.48367);
}

double NNtest::synapse0x25ec230() {
   return (neuron0x25ebb90()*42.9957);
}

