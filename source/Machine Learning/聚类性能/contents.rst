聚类性能
================


==========================
Ⅰ. Davies-Bouldin Index
==========================

It estimates the cohesion based on the distance from the points in a cluster to its centroid and the separation based on the distance between centroids. It is defined as:

:math:`{\rm DBI}(\mathbb{C})=\dfrac{1}{K}\sum\limits_{C_k\in\mathbb{C}}\max\limits_{C_j\in\mathbb{C}\backslash C_k}\left\{ \dfrac{\dfrac{1}{|C_k|}\sum\limits_{x_i\in C_k}d(x_i,\overline{C_k})+\dfrac{1}{|C_j|}\sum\limits_{x_i\in C_j}d(x_i,\overline{C_j})}{d(\overline{C_k},\overline{C_j})} \right\}`


==========================
Ⅱ. Dunn Index
==========================

It is a ratio-type index where the cohesion is estimated by the nearest neighbor distance and the separation by the maximum cluster diameter. The original index is defined as:

:math:`{\rm DI}(\mathbb{C})=\dfrac{\min\limits_{C_k\in\mathbb{C}}\left\{ \min\limits_{C_j\in\mathbb{C}\backslash C_k}\left[ \min\limits_{x_i\in C_k}\min\limits_{x_j\in C_j}\{ d(x_i,x_j) \} \right] \right\}}{\max\limits_{C_k\in\mathbb{C}}\left\{ \max\limits_{x_i,x_j\in C_k}\left[d(x_i,x_j)\right] \right\} }`

:math:`{\rm where}\ \ \min\limits_{x_i\in C_k}\min\limits_{x_j\in C_j}\{ d(x_i,x_j)\}=\min\limits_{x_i\in C_k;x_j\in C_j}\{ d(x_i,x_j)\}`



==========================
Ⅲ.  Silhouette Index
==========================

This index is a normalized summation-type index. The cohesion is measured based on the distance between all the points in the same cluster and the separation is based on the nearest neighbor distance. It is defined as:

:math:`{\rm Sil}(\mathbb{C})=\dfrac{1}{N}\sum\limits_{C_k\in\mathbb{C}}\sum\limits_{x_i\in C_k}\dfrac{\min\limits_{C_j\in\mathbb{C}\backslash C_k}\left[ \dfrac{1}{|C_j|}\sum\limits_{x_j\in C_j}d(x_i,x_j) \right]-\dfrac{1}{|C_k|}\sum\limits_{x_j\in C_k}d(x_i,x_j)}{\max\left\{ \dfrac{1}{|C_k|}\sum\limits_{x_j\in C_k}d(x_i,x_j),\min\limits_{C_j\in\mathbb{C}\backslash C_k}\left[ \dfrac{1}{|C_j|}\sum\limits_{x_j\in C_j}d(x_i,x_j) \right] \right\}}`

:math:`{\rm where}\ \ \sum\limits_{C_k\in\mathbb{C}}\sum\limits_{x_i\in C_k}\Longleftrightarrow \sum\limits_{C_k\in\mathbb{C};x_i\in C_k}`