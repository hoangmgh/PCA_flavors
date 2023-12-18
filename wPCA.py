from scipy.sparse.linalg import svds
from scipy.sparse import csr_matrix
import numpy as np
import scanpy as sc
def weighted_pca(anndata,weights,n_comps, corr=True,
                 max_value=10):
   '''
        This function calculate a weighted version of PCA as suggested in https://doi.org/10.1016/j.medj.2022.05.002
        by Korsunsky et al
        input:
                anndata: a log10CPM normalized count matrix. This should not be scaled in any ways
                n_comps: number of components for wPCA

    '''
    from scipy.sparse.linalg import svds
    from scipy.sparse import csr_matrix
    import numpy as np
    assert anndata.shape[0] == weights.shape[0]
    if max_value is not None:
        assert max_value > 0
    def vars(a, axis=None):
        """ Variance of sparse matrix a
        var = mean(a**2) - mean(a)**2
        """
        a_squared = a.copy()
        a_squared.data **= 2
        return a_squared.mean(axis) - np.square(a.mean(axis))

    def stds(a, axis=None):
        """ Standard deviation of sparse matrix a
        std = sqrt(var(a))
        """
        return np.sqrt(vars(a, axis))

    wanndataX=anndata.X.copy()

    sklearn.utils.sparsefuncs.inplace_row_scale(wanndataX, weights)
    mu=wanndataX.mean(axis=0)
    sig=np.array((1/stds(wanndataX,axis=0)).T)
    ##### We want to scale the old data (not the new data using mean and std of the new data):
    print("shifting by weighted mean:")
    wanndataX=scipy.sparse.csr_matrix(anndata.X.copy()-mu)
    print("row scaling with weighted standard deviation")
    sklearn.utils.sparsefuncs.inplace_column_scale(wanndataX, sig)

    del mu,sig
    if corr=True:
        #####And then run scaling on this matrix one more time observation-wise:
        mu=wanndataX.mean(axis=1)
        print(mu.shape)
        sig=np.array((1/stds(wanndataX,axis=1)))
        print(sig.shape)
        wanndataX=scipy.sparse.csr_matrix(wanndataX-mu)
        print(" scaling with weighted standard deviation")

        sklearn.utils.sparsefuncs.inplace_row_scale(wanndataX, sig)
        if max_value is not None:
            wanndataX[wanndataX>=max_value]=max_value
            wanndataX[wanndataX<=-1*max_value]=-1*max_value

    ##### and then multiply this column sqrt of weights:
    print("adding weights:")
    sklearn.utils.sparsefuncs.inplace_row_scale(wanndataX, np.sqrt(weights))
    del mu,sig
    print("running SVDs:")
    u, s, v = svds( wanndataX.T, k=n_comps)
    v=scipy.sparse.csr_matrix(v)
    


    ############ 
    print("outputing V and U:")
    w=np.array((1/np.sqrt(weights)).T)
    sklearn.utils.sparsefuncs.inplace_column_scale(v,w)
    sklearn.utils.sparsefuncs.inplace_row_scale(v,np.array(s))

    
    gc.collect()
    return np.array(v.todense().T),np.array(u)

def generate_weights(anndata,Batch_key):
    assert Batch_key in anndata.obs.columns
    freq= anndata.obs[Batch_key].value_counts().reset_index().rename({"index":Batch_key,Batch_key:"count"},axis=1)
    freq["count"]=freq["count"]/np.sum(freq["count"])
    freq=pd.merge(dt.obs[[Batch_key]],freq,how="left",on=Batch_key)
    w=1/freq["count"]
    w=w/np.sum(w)
    return(w)
