import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans


def data_generator(n ,dim,noise_v, number,plot_d = True):
    K = len(n)
    center = np.array([])
    dim_zero = np.random.choice(np.arange(dim),noise_v, replace=False)
    for i in range(dim):
        coin = np.random.rand(K) < 0.5
        col = np.where(coin,np.random.uniform(-6, -3, K),np.random.uniform(3, 6, K))
        center = np.append(center, np.array(col))
    center = center.reshape(K,-1)

    center[:,dim_zero] = 0
    Sigmar = np.eye(dim)
    all_data_normal = np.array([])
    for i in range(K):
        all_data_normal = np.append(all_data_normal,
                    np.random.multivariate_normal(center[i,:], Sigmar,n[i]))

    all_data_normal = all_data_normal.reshape(-1,dim)

    outliner_idx = np.random.choice(all_data_normal.shape[0],int(number*all_data_normal.shape[0]))

    r_1 = np.random.uniform(5,10,100)
    r_2 = np.random.uniform(-10,-5,100)
    r = np.append(r_1,r_2)
    gap_matrix = np.zeros(all_data_normal.shape[0]*dim).reshape(all_data_normal.shape[0],dim)
    for i in outliner_idx:
        gap_matrix[i,:] = np.random.choice(r,dim)
    oulier_dataset_all = all_data_normal + gap_matrix
    if plot_d == True:
        plt.scatter(oulier_dataset_all[:, 0], oulier_dataset_all[:,1],
                    c="navy", marker='o', label='Data without outlier')
        plt.show()
    return {"data":oulier_dataset_all,"outlier_idx":outliner_idx,"center":center}




#algorithm of ARSKC
def distance_point(x,y):
    return np.sqrt(((x - y)**2).sum())


def soft_thresholding_E(x,m,lam_1):
    return (x - m)* np.maximum(0, 1-(lam_1*(distance_point(x,m)**(-1))))

def soft_thresholding_v(x,lam_2):
    if x > lam_2:
        return x - lam_2
    if np.abs(x) <= lam_2:
        return 0
    if x < -lam_2:
        return x + lam_2


def W_rocover(x):
    if x == 0:
        return 1
    else:
        return x


def robust_kmean(dataset,lam_1,K,W):
    total_mean = dataset.mean(axis=0)
    all_distence = np.array([])
    dataset_n = dataset.shape[0]
    dataset_p= dataset.shape[1]
    for i in range(dataset_n):
        all_distence = np.append(all_distence,distance_point(dataset[i],total_mean))

    all_oulier_idx = np.argsort(all_distence)[(dataset_n - 5):dataset_n]

    innit_E = np.zeros((dataset_n,dataset_p))
    innit_E[all_oulier_idx,:] = dataset[all_oulier_idx,:]

    clu_dataset = dataset - innit_E
    weight_dataset = clu_dataset*W
    kmeans = KMeans(n_clusters=K)
    kmeans.fit(weight_dataset)
    old_center = kmeans.cluster_centers_
    old_cluster = kmeans.labels_



    iter = 0
    while True:
        lable_dict = {}
        for i in np.unique(old_cluster):
            lable_dict[i] = np.where(old_cluster == i)[0].tolist()


        E_n = innit_E
        for i in lable_dict:
            for j in lable_dict[i]:
                E_n[j,:] = soft_thresholding_E(dataset[j,:],old_center[i],lam_1 = lam_1)

        Kmeans = KMeans(n_clusters=K)
        Kmeans.fit(clu_dataset)

        new_center = Kmeans.cluster_centers_

        if distance_point(new_center,old_center) <= 1e-12 or iter == 200:
            break

        old_center = new_center

        old_cluster = Kmeans.labels_

        iter = iter + 1

    outlier_idex_output = np.where(np.diag(E_n@E_n.T) != 0)[0]

    return {
        "center": old_center,
        "label_dict": old_cluster,
        "E_est":E_n,
        "iter": iter,
        "outlier_index" :outlier_idex_output
    }



def find_W(W, dataset, E_res, cluster_result, cluster_center, lam_2):
    recover_W = np.array([W_rocover(i) for i in W])

    E_return = E_res / np.sqrt(recover_W)
    new_dataset = dataset - E_return


    all_mean = np.mean(new_dataset, axis=0)
    all_diff = np.sum((new_dataset - all_mean) ** 2, axis=0)

    lable_dict = {i: np.where(cluster_result == i)[0].tolist()
                  for i in np.unique(cluster_result)}

    cluster_diff = []
    for i in lable_dict:
        cluster_subset = new_dataset[lable_dict[i], :]
        cluster_center_i = np.mean(cluster_subset, axis=0)
        cluster_diff.append(np.sum((cluster_subset - cluster_center_i) ** 2, axis=0).tolist())

    cluster_diff_K = np.sum(cluster_diff, axis=0)
    Q_j = all_diff - cluster_diff_K
    after_thresholding_Q_j = np.array([soft_thresholding_v(j, lam_2=lam_2) for j in Q_j])
    W_n = after_thresholding_Q_j / np.sqrt(np.sum(after_thresholding_Q_j ** 2))

    return W_n



def total_ARSKC(dataset, lam_1, lam_2, K):
    w_init = np.ones(dataset.shape[1]) / np.sqrt(dataset.shape[1])
    W = w_init
    while True:
        one_result = robust_kmean(dataset=dataset, W=W, lam_1=lam_1, K=K)

        W_n_1 = find_W(W=W,dataset=dataset,E_res=one_result["E_est"],cluster_result=one_result["label_dict"],
                       cluster_center=one_result["center"],lam_2=lam_2)

        judge_w = np.sum(np.abs(W_n_1 - W)) / np.sum(np.abs(W))

        print(iter_W)

        if judge_w <= 1e-1 or iter_W == 50:
            break

        W = W_n_1

    return {"weight":W_n_1, "center":one_result["center"],"label":one_result["label_dict"]}



# data_g_1 = data_generator([10,10],dim=10,noise_v=2,number=0.1,plot_d=False)
# dataset = data_g_1["data"]
# total_ARSKC(dataset=dataset,lam_1 = 20,lam_2 = 1.1,K = 2)


###############gpp generated by LLM maynot correct

def GetWCSS(x, Cs, ws=None):
    wcss_perfeature = np.zeros(x.shape[1])
    
    for k in np.unique(Cs):
        whichers = (Cs == k)
        if np.sum(whichers) > 1:
            cluster_data = x[whichers, :]
            centered_data = cluster_data - np.mean(cluster_data, axis=0)
            wcss_perfeature += np.sum(centered_data**2, axis=0)
    centered_x = x - np.mean(x, axis=0)
    total_ss_perfeature = np.sum(centered_x**2, axis=0)
    bcss_perfeature = total_ss_perfeature - wcss_perfeature
    
    result = {
        'wcss_perfeature': wcss_perfeature,
        'wcss': np.sum(wcss_perfeature),
        'bcss_perfeature': bcss_perfeature
    }
    
    if ws is not None:
        result['wcss_ws'] = np.sum(wcss_perfeature * ws)
    
    return result


def cal_Gap(dataset, k, lambda_1, lambda_2):
    
    sample_run = total_ARSKC(dataset=dataset, lam_1=lambda_1, lam_2=lambda_2, K=k)
    
    if 'E_est' in sample_run:
        x_star = dataset - sample_run['E_est']
    else:
        x_star = dataset
    
    wcss_result = GetWCSS(x_star, sample_run['label'])
    o_all = np.sum(sample_run['weight'] * wcss_result['bcss_perfeature'])
    nperms = 25
    permx = []
    
    np.random.seed(42) 
    for i in range(nperms):
        permuted = np.zeros_like(dataset)
        for j in range(dataset.shape[1]):
            permuted[:, j] = np.random.permutation(dataset[:, j])
        permx.append(permuted)
    
    permtots = []
    for K in range(nperms):
        perm_out = total_ARSKC(dataset=permx[K], lam_1=lambda_1, lam_2=lambda_2, K=k)
        if 'E_est' in sample_run:
            perm_data = dataset - sample_run['E_est']
        else:
            perm_data = dataset
            
        perm_bcss = GetWCSS(perm_data, sample_run['label'])['bcss_perfeature']
        permtots.append(np.sum(perm_out['weight'] * perm_bcss))
    
    t_iter = sample_run.get('t_iter', 0)
    okm_it = sample_run.get('okm_it', 0)
    
    if (t_iter >= 15) or (okm_it >= 50):
        Gap = None  # Python中用None代替R的NA
    else:
        permtots = np.array(permtots)
        if np.any(permtots <= 0):
            print("Warning: Some permutation totals are <= 0, Gap calculation may be invalid")
            Gap = None
        else:
            Gap = np.log(o_all) - np.mean(np.log(permtots))
    
    return {
        'Gap': Gap,
        'k': k,
        'lambda_1': lambda_1,
        'lambda_2': lambda_2,
        'B_model': sample_run
    }













