## Introduction
This is the MATLAB code for "<a href="https://arxiv.org/abs/2112.13199">A Spectral Method for Joint Community Detection and Orthogonal Group Synchronization</a>" by Yifeng Fan, Yuehaw Khoo, and Zhizhen Zhao. 

## Getting Started
In MATLAB, Please run 
~~~
demo.m 
~~~
which contains an example of recovering the cluster memberships and orthogonal transformations based on our model. Feel free to adjust the parameters such as the number of clusters and cluster sizes.

## Blockwise column pivoted QR factorization
An important ingredient of our algorithm is the **blockwise column pivoted QR factorization** (blockwise CPQR) proposed in this project, which might be also useful on other applications. If necessary, please find the implementation in
~~~
functions/qr_householder_block.m
~~~

## Miscellaneous
If you find this repository useful for your research, please consider citing our paper:

    @article{fan2021spectral,
      title={A Spectral Method for Joint Community Detection and Orthogonal Group Synchronization},
      author={Fan, Yifeng and Khoo, Yuehaw and Zhao, Zhizhen},
      journal={arXiv preprint arXiv:2112.13199},
      year={2021}
    }

## Miscellaneous

Please send any problems/questions you might have to <yifengf2@illinois.edu>.
