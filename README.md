# ResFouriONet
Residual Fourier Operator Networks for 2D transient bioheat transfer problems with parametric laser excitation

## Abstract

Accurate modeling of laser-tissue interactions is essential for applications such as multiphoton brain imaging and laser-based therapeutics, where precise thermal control is critical. However, traditional numerical solvers for bioheat transfer, especially when paired with Monte Carlo (MC) simulations for light transport, are computationally prohibitive for real-time use. While operator learning offers a scalable alternative, existing architectures often struggle with extrapolatory tasks and rely on large, labeled datasets that are costly to generate. In this study, we introduce ResFouriONet, a residual deep operator network tailored for modeling 2D transient bioheat transfer with parametric laser sources. Our architecture features a dual-pathway branch network comprising a Fourier sub-branch for capturing global thermal patterns and a cascaded residual sub-branch for localized heating, fused via a self-attention mechanism. To reduce the reliance on expensive MC data, we introduce a synthetic laser source function generation scheme that combines parameterized Gaussian profiles, exponential decay, and random field perturbations to emulate the complexity of MC-derived fluence maps. ResFouriONet achieves a 5× reduction in prediction error compared to vanilla deep operator networks, despite having fewer parameters. It attains an average relative error of 0.93% on unseen synthetic data and 1.45% on 630 zero-shot MC-generated source functions. Together, the proposed architecture and data generation strategy enable accurate, real-time prediction of laser-induced heating, paving the way for simulation-free planning and control in biomedical imaging and laser-based therapies.

## ResFouriONet architecture

We have drawn inspiration from the works of DeepONet, FNO, ResNet, U-Net, and Network in Network in designing the network architecture of ResFouriONet. 

<img width="2480" height="2185" alt="architecture_github" src="https://github.com/user-attachments/assets/e682e01b-55bb-43da-8b73-54925cd8df40" />

## Generalization performance 

A comparison of prediction performance between ResFouriONet, vanilla DeepONet and FNO-3D for two representative test cases is shown in panels (a) and (b). Each panel displays the ground-truth temperature field at the final time point (left), the corresponding network prediction (middle), and the absolute error (right). ResFouriONet captures both global thermal distributions and localized heating features with high fidelity, yielding negligible errors outside high-gradient regions and low peak-temperature deviations (average error = 1.37% across 200 test cases).

<img width="2244" height="3366" alt="Fig3_v2_github" src="https://github.com/user-attachments/assets/201598d5-157c-4735-a816-096a862ab88f" />

## Zero-shot testing performance

We evaluated the zero-shot performance of ResFouriONet on two representative MC-derived source functions. The parameter set used for (a) was FOV = 0.1 mm, NA = 0.95, imaging depth = 6l_t, average power = 500 mW, wavelength = 920 nm, and that for (b) was FOV = 0.5 mm, NA = 1.05, imaging depth = 2l_t, average power = 200 mW, wavelength = 1,035 nm. In both cases, ResFouriONet accurately reproduced peak temperatures and global thermal distributions compared to finite-difference ground truth solutions derived from MC fluence inputs (see MATLAB code provided). Across the full test set of 630 MC-generated functions, the network achieved an average testing error of 1.45%, demonstrating strong extrapolatory capability despite being trained solely on synthetic laser source functions.

<img width="2480" height="1299" alt="Fig4_draft_v2_github" src="https://github.com/user-attachments/assets/c4839d6f-5908-4505-afe0-da2be189226e" />

\n \n

The following box and violin plots show the zero-shot inference performance of ResFouriONet across 630 MC-derived source functions, illustrating the dependence of predictive accuracy on key optical and geometrical parameters. 

<img width="2480" height="1512" alt="Fig5_draft_v2_box_github" src="https://github.com/user-attachments/assets/7cc3eac3-969c-429f-be98-6b92d1deab72" />

