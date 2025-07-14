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



data_g_1 = data_generator([10,10],dim=10,noise_v=2,number=0.1,plot_d=False)
dataset = data_g_1["data"]



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



w_init = np.ones(dataset.shape[1]) / np.sqrt(dataset.shape[1])
result_sample_run_1 = robust_kmean(dataset = dataset, lam_1 = 20,K = 2, W = w_init)




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



def total_ARSKC(dataset, lam_1, lam_2, K, max_iter=400, tolerance=1e-1):
    w_init = np.ones(dataset.shape[1]) / np.sqrt(dataset.shape[1])
    W = w_init
    for iter_W in range(max_iter):
        one_result = robust_kmean(dataset=dataset, W=W, lam_1=lam_1, K=K)

        W_n_1 = find_W(
            W=W,
            dataset=dataset,
            E_res=one_result["E_est"],
            cluster_result=one_result["label_dict"],
            cluster_center=one_result["center"],
            lam_2=lam_2
        )

        judge_w = np.sum(np.abs(W_n_1 - W)) / np.sum(np.abs(W))
        print(iter_W)

        if judge_w < tolerance:
            break

        W = W_n_1

    return {"weight":W_n_1, "center":one_result["center"],"label":one_result["label_dict"]}



# total_ARSKC(dataset=dataset,lam_1 = 20,lam_2 = 2,K = 2)









