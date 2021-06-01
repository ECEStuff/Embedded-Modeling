# Embedded-Modeling
Embedding System Modeling using the Canny Edge Detector. Modeled in SpecC.

Versions
- canny.c : The original Canny Edge Detector by Mike Heath, with a bug fix of fixing a off-by-1 error in non_max_supp.
- canny_testonly.cpp : Performance estimation of the Canny Edge Detector.
- canny_v2.sc : Initial version in SpecC SLDL. Replaced dynamic allocation and user-adjustable configuration parameters with static allocation and hard-coded parameters.
- canny_v3.sc : Added test bench. Also converted from single image processing to real-time video processing.
- canny_v4.sc : Modularized the DUT. Previously, DUT contained local functions; these functions are now standalone modules.
- canny_v5.sc : Parallelization and pipelining of the DUT: pipelined the modules inside the DUT and parallelized BlurX and BlurY.
- canny_v6.sc : Optimization using different compiler optimization (O2, O3) and simulated hardware (replacing measurements on a simulated Raspberry Pi 3 with Raspberry Pi 4). 
