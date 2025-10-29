# ResFouriONet
Residual Fourier Operator Networks for 2D transient bioheat transfer problems with parametric laser excitation

## Abstract

Accurate modeling of laser-tissue interactions is essential for applications such as multiphoton brain imaging and laser-based therapeutics, where precise thermal control is critical. However, traditional numerical solvers for bioheat transfer, especially when paired with Monte Carlo (MC) simulations for light transport, are computationally prohibitive for real-time use. While operator learning offers a scalable alternative, existing architectures often struggle with extrapolatory tasks and rely on large, labeled datasets that are costly to generate. In this study, we introduce ResFouriONet, a residual deep operator network tailored for modeling 2D transient bioheat transfer with parametric laser sources. Our architecture features a dual-pathway branch network comprising a Fourier sub-branch for capturing global thermal patterns and a cascaded residual sub-branch for localized heating, fused via a self-attention mechanism. To reduce the reliance on expensive MC data, we introduce a synthetic laser source function generation scheme that combines parameterized Gaussian profiles, exponential decay, and random field perturbations to emulate the complexity of MC-derived fluence maps. ResFouriONet achieves a 5× reduction in prediction error compared to vanilla deep operator networks, despite having fewer parameters. It attains an average relative error of 0.93% on unseen synthetic data and 1.45% on 630 zero-shot MC-generated source functions. Together, the proposed architecture and data generation strategy enable accurate, real-time prediction of laser-induced heating, paving the way for simulation-free planning and control in biomedical imaging and laser-based therapies.

## ResFouriONet architecture

We have drawn inspiration from the works of DeepONet, FNO, ResNet, U-Net, and Network in Network in designing the network architecture of ResFouriONet. 

<img width="468" height="412" alt="image" src="https://github.com/user-attachments/assets/647b8b50-e764-44d5-8075-dcb893b96d64" />
