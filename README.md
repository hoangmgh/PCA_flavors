
In [Cross-tissue, single-cell stromal atlas identifies shared pathological fibroblast phenotypes in four chronic inflammatory diseases](https://doi.org/10.1016/j.medj.2022.05.002), 
the authors use a technique called "weighted PCA", together with harmony, 
to remove batch effect across a wide variety of diseases. 
The most important observation they made is the stark difference between the number of cells between tissues. 
To dive into this, I would like to explain the concept of weighted PCA, first by formalizing  the idea of weighted expectation, weighted variance, and eventually weighted covariance matrix. While there are various implementation of weighted PCA out there, this is perhaps the easiest implementation, most intuitive, and also well generalized from the original definition of PCA. 
 First, we can define a weighted inner-product in the Euclidean space: 

$$
    \langle x,y \rangle_W = x^T W y  
$$

where the diagonal entries of $W$ stores the weights ($diag(W)=\vec{w}$) and their entries must sum to $1$. You can check that this is indeed an inner product by checking its properties. As a result, the weighted norm simply follows:

$$
    \|x\|_W = \langle x,x \rangle_W = (x^T W x)^{1/2} 
$$
Then, we  can define a weighted mean of a vector $ \vec{x} \in \mathbb{R^n} $. You can think of a mean as a dot-product as well!

$$
    \mu_x^W  = \vec{x}^T \cdot \vec{w}  = \vec{x}^T \cdot W \cdot  \vec{1} = \langle x, \vec{1} \rangle_W 
$$

and $diag(W)=w$. You can always replace the weight vector with the diagonal matrix $W$ of the same size!
In the unweighted case, we simply have all entries of $w$ to be $\frac{1}{n}$.
Now, we can define weighted covariance  of $x$ and $y$ as:

\begin{align}
    Cov(x,y)_W = \langle x-\mu_x^W ,y-\mu_y^W \rangle_W  = (x-\mu_x^W)^T W (y-\mu_y^W)
\end{align}
In short, most of our measures, i.e. correlation, covariance, mean, variance, are replaced with the weighted version. I think it makes sense that this has to be built from the ground up using a different version of the dot-product.
The weights for each observation can also be interpreted as corresponding to the frequency of each observation. In an imbalance situation, it is favorable to  incorporate this weight to reflect the frequency of different classes of observation. 

Now, for a large matrix $A$ of form $\mathbb{R^{g \times c }}$, we can "bulk" compute the sample covariance matrix the following way. First, 
center each gene at the weighted mean and inversely scale them by the weighted standard deviation. From there, we can scale the matrix observation-wise
 (so that later on $AA^T$  actually sample correlation matrix):

The covariance matrix is, in fact, no longer $AA^T$ but $AWA^T$, due to our definition of the covariance above. The weighted PCA from here can be rewritten as diagonalizing (eigendecomposition):
$$
AWA^T = AW^{1/2} W^{1/2}A^T = AW^{1/2} (AW^{1/2})^T
$$
Here $W$ is diagonal so $W^{1/2}$ is the same as its transpose.
Therefore, diagonalizing $AWA^T$ is equivalent  to running SVD for  $AW^{1/2}$. We can then write $AW^{1/2}$ as
$$
AW^{1/2}=USV
$$
and hence
$$
A = USVW^{1/2}
$$
Under the new orthogonal basis spanned by $U$, the coordinates are   now given by  $SVW^{1/2}$ 
