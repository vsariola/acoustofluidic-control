# acoustofluidic-control

This repository contains code for TUNI acoustofluidic manipulation
setup. The setup is documented in the following papers:

Kyriacos Yiannacou and Veikko Sariola, "Controlled Manipulation and
Active Sorting of Particles Inside Microfluidic Chips Using Bulk
Acoustic Waves and Machine Learning", Langmuir, 37 (14), 4192–4199,
2021.

Kyriacos Yiannacou, Vipul Sharma, and Veikko Sariola, "Programmable
Droplet Microfluidics Based on Machine Learning and Acoustic
Manipulation", Langmuir, 38 (38), 11557–11564, 2022.

Kyriacos Yiannacou and Veikko Sariola, "Acoustic Manipulation of
Particles in Microfluidic Chips with an Adaptive Controller that Models
Acoustic Fields", Advanced Intelligent Systems, 5 (9), 2300058, 2023.

The code was written and tested in Matlab 2019a&b and Matlab 2021a.

Probably the most interesting part of the repository is the script [src/bandit_ctrl.m](src/bandit_ctrl.m),
which contains the implementation of our e-greedy and UCB1 control
algorithms, and [src/sfb_ctrl.m](src/sfb_ctrl.m), which contains the
script for AMA control algorithm (originally was about to be called
"sparse fourier basis controller" but got renamed to AMA controller
during the writing of the paper).

Also [src/+experiments/](src/+experiments/) is of interest, as it
contains basic scripts to run the single-particle, multi-particle
manipulation experiments, scripts to learn the sparse regression
controller, and a script to write the university logo.

## How to run tests

To test the code, run `runtests tests`.

## License

[MIT](LICENSE)