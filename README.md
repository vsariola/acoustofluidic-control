# acoustofluidic-control

This repository contains code for TUNI acoustofluidic manipulation
setup. The setup is documented in the following paper:

Kyriacos Yiannacou and Veikko Sariola, "Controlled Manipulation and
Active Sorting of Particles Inside Microfluidic Chips Using Bulk
Acoustic Waves and Machine Learning", Langmuir, 37 (14), 4192–4199,
2021.

The code was written and tested in Matlab 2019a&b and Matlab 2021a.

Probably the most interesting part of the repository is the script
[src/bandit_ctrl.m](src/bandit_ctrl.m), which contains the
implementation of our e-greedy and UCB1 control algorithms, and
[src/sfb_ctrl.m](src/sfb_ctrl.m), which contains the script for sparse
regression control algorithm.

Also [src/+experiments/](src/+experiments/) is of interest, as it
contains basic scripts to run the single-particle, multi-particle
manipulation experiments, scripts to learn the sparse regression
controller, and a script to write the university logo.

## How to run tests

To test the code, run `runtests tests`.

## License

[MIT](LICENSE)