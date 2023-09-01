# SIR-diffusion-in-two-layer-network

ganger_demo.m is an example that can convert multiple sets of time series data into a pairwise ganger causality matrix.
Epidemic_simulation_SIR original.m can realize SIR diffusion simulation on a single-layer network. Input: relationship matrix, infection coefficient, recovery rate. Output: SIR infection curve.
Epidemic_simulation_SIR_two_layer.m can implement the SIR diffusion model on a two-layer network. Input: two relationship matrices, coupling nodes, two infection coefficients, two recovery rates. Output: SIR infection curve corresponding to the two network layers.

Users can prepare the relationship matrix by themselves, or use the ganger relationship matrix in this directory, or find the generation functions of SW and BW random networks in Epidemic_simulation_SIR_two_layer.

This program is a simulation experiment of the paper. If you find it helpful, please cite the following two articles:
[1]Huang Q A, Zhao J C, Wu X Q. Financial risk propagation between Chinese and American stock markets based on multilayer networks[J]. Physica A: Statistical Mechanics and its Applications, 2022, 586: 126445.
[2]Jun-Chan Z, Qi-An H, Xiao-Qun W, et al. The impact of trade war on Shanghai stock exchange industry based on Granger causality network[J]. ACTA PHYSICA SINICA, 2021, 70(7) .

# 双层网络上的SIR传染病扩散

ganger_demo.m是一个能将多组时间序列数据转化为两两ganger因果关系矩阵的示例。
Epidemic_simulation_SIR original.m能实现单层网络上的SIR扩散模拟。输入：关系矩阵，感染系数，恢复率。输出：SIR感染曲线。
Epidemic_simulation_SIR_two_layer.m能实现双层网络上的SIR扩散模型。输入：两个关系矩阵，耦合节点，两个感染系数，两个恢复率。输出：两个网络层对应的SIR感染曲线。

使用者可以自行准备关系矩阵，也可以使用本目录下的ganger关系矩阵，也可以在Epidemic_simulation_SIR_two_layer里面找到SW和BW随机网络的生成函数。

本程序为论文的模拟实验。如果你觉得对你有帮助请引用下面两篇文章:
[1]Huang Q A, Zhao J C, Wu X Q. Financial risk propagation between Chinese and American stock markets based on multilayer networks[J]. Physica A: Statistical Mechanics and its Applications, 2022, 586: 126445.
[2]Jun-Chan Z, Qi-An H, Xiao-Qun W, et al. The impact of trade war on Shanghai stock exchange industry based on Granger causality network[J]. ACTA PHYSICA SINICA, 2021, 70(7).
