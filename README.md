# Non negative Matrix Factorization with prior and marking

In the context of high throughput RNA sequencing, we aim at solving the inverse problem of subtype mixing within a sample.

We formulate this problem as a matrix factorization problem, where the data matrix is factorized into two matrices, one representing the subtypes and the other the mixing proportions.

Let $X \in \mathbb{R}^{m \times n}$ be a count matrix where $n$ is the number of samples and $m$ is the number of genes. We aim at determining $K \in \mathbb{N}^*$ prototypes.

We assume that we have a matrix $M \in \{0,1\}^{m \times k}$ where $k\leq K$ is the number of prototypes on which we put a prior. That is, an entry is marked as 1 if we want this gene to represent, at least partially, the corresponding prototype. We aim at solving the following problem:

$ W^*,H^* = \arg\min_{H \in \mathbb{R}^{n \times K}_{+} ; W \in \mathbb{R}^{m \times K}} \text{dst}_X(WH) + \lambda g_M(W)$

where we have:

$$
\text{dst}_X(WH) = \frac{1}{2mn}|| WH - X||_F^2
$$

$$
g_M(W) =  \frac{1}{2mk} ||W[:,:k] \odot (1-M)||_F^2 = \frac{1}{2mk}\sum_{i=1}^{m} \sum_{j=1}^{k} W_{ij}^2 (1-M_{ij}) 
$$

To do so we use alterned gradient descent.

For more details and experiments, please refer to the [notebook](https://github.com/Theaublanchard/NMFMarkedtree/master/markedNMF.ipynb).
