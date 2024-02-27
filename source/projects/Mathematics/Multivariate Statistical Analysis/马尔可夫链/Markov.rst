马尔可夫链
====================

|


Ⅰ. 背景 --> Ⅱ. 原理 --> Ⅲ. 示例 --> Ⅳ. Python 代码

|



================
Ⅰ. 背景
================

“The future is independent of the past given the present” --- Andrey Andreyevich Markov.

马尔可夫链（Markov Chain）是机器学习和人工智能的基石. 马尔可夫链的核心思想是，过去所有的信息都蕴含在当下，未来只基于此时此刻而和过去无关. 虽然在哲学角度，这种极端的说法不具有辩证性，但马尔可夫链在很多时间序列模型中得到广泛的应用.



|

================
Ⅱ. 原理
================

定义序列状态为 :math:`X^{(i)},\ i=0,1,2,\ldots,n`. 例如，第 0 天是晴天，第 1 天是雨天，第 2 天是雪天.

满足马尔可夫假设：

.. math:: \Pr(X^{(n+1)}\mid X^{(n)},X^{(n-1)},\ldots,X^{(0)})=\Pr(X^{(n+1)}\mid X^{(n)})

第 n+1 次出现某状态的概率独立于先前所有状态，等价于在第 n 次状态发生后第 n+1 次的概率. 这种现象称之为无后效性、无记忆性.

统计以前所有相同状态转变的频率，得到状态转移矩阵（经济上也称支付矩阵），代表状态转变发生的可能性：

.. math::
    \begin{align}
        \bm{P}=
        \left[
        \begin{matrix}
            p_{11} & p_{12} & \cdots & p_{1r} \\
            p_{21} & p_{22} & \cdots & p_{2r} \\
            \vdots & \vdots & \ddots & \vdots \\
            p_{r1} & p_{r2} & \cdots & p_{rr} \\
        \end{matrix}
        \right]
        \notag
    \end{align}

其中，:math:`p_{ij}` 代表状态 :math:`X^{(i)}` 转变到状态 :math:`X^{(j)}` 发生的频率.

如果此时此刻各状态发生的概率所构成的向量已知：

.. math:: \bm{\vec{\pi}}(n)=[\pi_1(n),\pi_2(n),\ldots,\pi_r(n)]

那么，接下来的状态发生概率向量则可以表示为：

.. math:: \bm{\vec{\pi}}(n+1)=\bm{\vec{\pi}}(n)\bm{P}

也可以递归表出为：

.. math:: \bm{\vec{\pi}}(n+1)=\bm{\vec{\pi}}(0)\bm{P}^{n+1}




|

================
Ⅲ. 示例
================

三位选手在比赛，裁判记录了得分情况如下：


.. list-table:: 
  :widths: 2 2 2 2 2 2 2 2 2 2
  :header-rows: 1

  * - 局号
    - 1
    - 2
    - 3
    - 4
    - 5
    - 6
    - 7
    - 8
    - 9
  * - 得分
    - 3
    - 3
    - 2
    - 3
    - 1
    - 2
    - 1
    - 3
    - 1


例如第 7 局是选手 1 得分. 那么，状态转移链就可以写成：

.. math:: 332312131

转移矩阵则统计为：

.. math::
    \begin{align}
        \bm{P}=
        \begin{matrix}
             & 1 & 2 & 3 \\
            \hline
            1\ \vline & 0 & 1 & 1 \\
            2\ \vline & 1 & 0 & 1 \\
            3\ \vline & 2 & 1 & 1 \\
        \end{matrix}
        \Longleftrightarrow
        \begin{matrix}
             & 1 & 2 & 3 \\
            \hline
            1\ \vline & 0 & 0.5 & 0.5 \\
            2\ \vline & 0.5 & 0 & 0.5 \\
            3\ \vline & 0.5 & 0.25 & 0.25 \\
        \end{matrix}
        \notag
    \end{align}

例如选手 3 转移到选手 1 一共发生 2 次，即有 2 次都是选手 3 输，随之选手 1 胜.

假设第十局三位选手得分情况未知，权且都标记为 33% 的概率得分，那么：

.. math:: \bm{\vec{\pi}}(11)=\left[\dfrac{1}{3},\ \dfrac{1}{3},\ \dfrac{1}{3}\right]\times\bm{P}=\left[\dfrac{1}{3},\ \dfrac{1}{4},\ \dfrac{5}{12}\right]

简言之，接下来选手 1 有 33% 的概率得分，选手 2 有 25% 的概率得分，选手 3 有 41.67% 的概率得分.



|

================
Ⅳ. 代码
================


.. code-block:: Python
    :caption: Markov.py
    :emphasize-lines: 2,9,10,92
    :linenos:

    '''
    # System --> Windows & Python3.10.0
    # File ----> Markov.py
    # Author --> Illusionna
    # Create --> 2024/02/02 22:45:54
    '''
    # -*- Encoding: UTF-8 -*-

    # numpy==1.26.3
    # pandas==2.0.3

    import os
    import numpy as np
    import pandas as pd
    from typing import Literal

    def cls() -> None:
        os.system('')
        os.system('cls')
    cls()


    class MARKOV:
        """
        马尔可夫类.
        """
        def __init__(self, statuslink:Literal['状态链字符串, 例如: abbabaacaa']) -> None:
            """
            初始化构造函数.
            """
            self.__statuslink = statuslink

        def __CountAdjacentStatusFrequency(
                statuslink:str,
                adjacentStatus:Literal['例如: 相邻状态字符串 ab, 表示 a -> b']
        ) -> int:
            """
            私有函数: 计算状态链中相邻状态 adjacentStatus 出现的频数.
            """
            adjacentStatusNumber = 0
            for i in range(0, len(statuslink)-1, 1):
                if statuslink[i:i+2] == adjacentStatus:
                    adjacentStatusNumber = adjacentStatusNumber + 1
            return adjacentStatusNumber

        def PayOffMatrix(self, probability:bool=False) -> pd.DataFrame:
            """
            公有函数: 计算支付矩阵.
            """
            uniqueItems = np.unique(list(self.__statuslink))
            df = pd.DataFrame(index=uniqueItems, columns=uniqueItems)
            for i in uniqueItems:
                for j in uniqueItems:
                    df.loc[i, j] = MARKOV.__CountAdjacentStatusFrequency(self.__statuslink, i+j)
            if probability == True:
                df = df.div(df.sum(axis=1), axis='index')
            self.__df = df.div(df.sum(axis=1), axis='index')
            return df

        @staticmethod
        def NextTransitionProbability(
                payOffMatrixProbability:pd.DataFrame,
                occurrenceProbability:dict|Literal['状态发生概率, 如字典形式 {a:0.5, b:0.3, c:0.2}'],
                step:int|Literal['接下来第 step 步, 默认为 1 步'] = 1,
                significant:int|Literal['小数点精度, 自然数, 如果异常, 可降低'] = 4
        ) -> np.array:
            """
            静态函数: 计算接下来的出现概率.
            """
            if len(occurrenceProbability.keys()) != payOffMatrixProbability.shape[0]:
                assert print(f'\033[3m\033[33m转移矩阵 {payOffMatrixProbability.shape[0]}x{payOffMatrixProbability.shape[0]} 维, 状态发生概率字典向量 1x{len(occurrenceProbability.keys())} 维, 无法相乘.\033[0m')
            tmp = 0
            L = []
            for (key, value) in occurrenceProbability.items():
                tmp = tmp + value
                L.append(value)
            try:
                np.testing.assert_approx_equal(tmp, 1.0, significant=significant)
            except:
                assert print(f'\033[3m\033[33m状态发生概率字典向量求和 {tmp} 不近似 1, 检查字典向量正确性或降低 significant={significant} 精度值\033[0m')
            del tmp
            output = np.array(L).dot(
                np.linalg.matrix_power(
                    np.array(payOffMatrixProbability, dtype=np.float64),
                    n = step
                )
            )
            return output


    if __name__ == '__main__':
        obj = MARKOV(statuslink='332312131')
        
        print(obj.PayOffMatrix(probability=False))

        nextTransitionProbability = MARKOV.NextTransitionProbability(
            payOffMatrixProbability = obj.PayOffMatrix(probability=True),
            occurrenceProbability = {
                '1': 1/3,
                '2': 1/3,
                '3': 1/3
            },
            step = 1    # 如果步长很大, 则收敛情况下, 趋于转移矩阵的极限分布.
        )
        print(nextTransitionProbability)