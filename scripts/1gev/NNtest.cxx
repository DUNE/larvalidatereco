#include "NNtest.h"
#include <cmath>

double NNtest::Value(int index,double in0,double in1,double in2) {
   input0 = (in0 - -1)/1;
   input1 = (in1 - -0.967455)/0.177442;
   input2 = (in2 - -0.904929)/0.293313;
   switch(index) {
     case 0:
         return neuron0x2a8cc50();
     default:
         return 0.;
   }
}

double NNtest::Value(int index, double* input) {
   input0 = (input[0] - -1)/1;
   input1 = (input[1] - -0.967455)/0.177442;
   input2 = (input[2] - -0.904929)/0.293313;
   switch(index) {
     case 0:
         return neuron0x2a8cc50();
     default:
         return 0.;
   }
}

double NNtest::neuron0x2a7aa70() {
   return input0;
}

double NNtest::neuron0x2a7adb0() {
   return input1;
}

double NNtest::neuron0x2a8bd40() {
   return input2;
}

double NNtest::input0x2a8c120() {
   double input = -0.243376;
   input += synapse0x2a7b180();
   input += synapse0x2a8c3d0();
   input += synapse0x2a8c410();
   return input;
}

double NNtest::neuron0x2a8c120() {
   double input = input0x2a8c120();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double NNtest::input0x2a8c450() {
   double input = -0.544843;
   input += synapse0x2a8c790();
   input += synapse0x2a8c7d0();
   input += synapse0x2a8c810();
   return input;
}

double NNtest::neuron0x2a8c450() {
   double input = input0x2a8c450();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double NNtest::input0x2a8c850() {
   double input = 0.118711;
   input += synapse0x2a8cb90();
   input += synapse0x2a8cbd0();
   input += synapse0x2a8cc10();
   return input;
}

double NNtest::neuron0x2a8c850() {
   double input = input0x2a8c850();
   return ((input < -709. ? 0. : (1/(1+exp(-input)))) * 1)+0;
}

double NNtest::input0x2a8cc50() {
   double input = 0.242197;
   input += synapse0x2a8ce70();
   input += synapse0x2a8ceb0();
   input += synapse0x2a8cef0();
   return input;
}

double NNtest::neuron0x2a8cc50() {
   double input = input0x2a8cc50();
   return (input * 1)+0;
}

double NNtest::synapse0x2a7b180() {
   return (neuron0x2a7aa70()*0.20476);
}

double NNtest::synapse0x2a8c3d0() {
   return (neuron0x2a7adb0()*-0.276663);
}

double NNtest::synapse0x2a8c410() {
   return (neuron0x2a8bd40()*-0.0622639);
}

double NNtest::synapse0x2a8c790() {
   return (neuron0x2a7aa70()*-0.0296275);
}

double NNtest::synapse0x2a8c7d0() {
   return (neuron0x2a7adb0()*-0.435347);
}

double NNtest::synapse0x2a8c810() {
   return (neuron0x2a8bd40()*0.176457);
}

double NNtest::synapse0x2a8cb90() {
   return (neuron0x2a7aa70()*-0.19548);
}

double NNtest::synapse0x2a8cbd0() {
   return (neuron0x2a7adb0()*0.227789);
}

double NNtest::synapse0x2a8cc10() {
   return (neuron0x2a8bd40()*0.314024);
}

double NNtest::synapse0x2a8ce70() {
   return (neuron0x2a8c120()*0.489482);
}

double NNtest::synapse0x2a8ceb0() {
   return (neuron0x2a8c450()*-0.221651);
}

double NNtest::synapse0x2a8cef0() {
   return (neuron0x2a8c850()*0.236219);
}

