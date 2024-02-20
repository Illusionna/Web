三次埃尔米特插值
================

.. role:: red
    :class: red

.. role:: blue
    :class: blue

.. role:: green
    :class: green

.. raw:: html

    <style>

        .red {
            color:red;
        }
        .blue {
            color:blue;
        }
        .green {
            color:green;
        }

    </style>

|

Ⅰ. 前言 --> Ⅱ. 作用 --> Ⅲ. 两种典型插值原理 --> Ⅳ. Python 代码


|

================
Ⅰ. 前言
================

涉及到部分拉格朗日插值基函数和牛顿插值均差系数，建议先阅读：


:doc:`../拉格朗日插值/Lagrange` ______ :doc:`../牛顿插值/Newton`



|

================
Ⅱ. 作用
================

埃尔米特（Hermite）插值是一类概念，可以认为不是某种具体的插值方法. 前面拉格朗日插值和牛顿插值只要求插值节点函数值与实际因变量相等即可，但在工程领域，像无人驾驶汽车的拐弯，这绝不是只要求函数值相等就可以了，如果曲率或者导数、导导数（二阶导数）甚至更高阶导数也能缓一点，那么转弯相对就平稳.

埃尔米特插值是要求插值函数与本源函数，在函数值角度相等，而且也在特点节点导数（或高阶导数）值也相等的一类概念. 为了避免龙格现象，工程上一般常用三次埃尔米特插值法已经能达到很好的效果.


|

=====================
Ⅲ. 两种典型插值原理
=====================

    1. :red:`三点三次`:green:`埃尔米特插值`

    自变量： :math:`(x_0 < x_1 < x_2)`

    因变量： :math:`(y_0,\ y_1,\ y_2)`

    满足条件：
        a. :math:`H(x_i)=y_i,\ \forall i\in\{0,1,2\}`

        b. :math:`H'(x_1)=y_1'`

    定义插值多项式【除三次项系数被修正，其余符合牛顿插值 :math:`f[x_0,x_1,\ldots,x_n]` 多项式形式】：

    .. math:: \displaystyle H(x)=f[x_0]\times 1+f[x_0,x_1](x-x_0)+f[x_0,x_1,x_2](x-x_0)(x-x_1)+C(x-x_0)(x-x_1)(x-x_2)

    由于在中间节点处导数相等，所以对 :math:`H(x)` 求导后代入 (b.) 条件可得：

    .. math:: C=\displaystyle\dfrac{y_1'-f[x_0,x_1]\times 1-f[x_0,x_1,x_2](x_1-x_0)}{(x_1-x_0)(x_1-x_2)}

    余项【其中 :math:`\xi=\xi(x)\in(x_0,x_2)`】：

    .. math:: R(x)=y(x)-H(x)=\dfrac{f^{(4)}(\xi)}{4!}(x-x_0)(x-x_1)(x-x_1)(x-x_2)


    2. :green:`两点三次`:blue:`埃尔米特插值`

    自变量： :math:`(x_0 < x_1)`

    因变量： :math:`(y_0,\ y_1)`

    满足条件：
        a. :math:`H(x_0)=y_0,\ H(x_1)=y_1`

        b. :math:`H'(x_0)=y_0',\ H'(x_1)=y_1'`

    定义插值多项式：

    .. math:: H(x)=y_0\alpha_0(x)+y_1\alpha_1(x)+y_0'\beta_0(x)+y_1'\beta_1(x)

    其中 :math:`\alpha(x)` 和 :math:`\beta(x)` 都是三次埃尔米特多项式基底, 且满足一定约束条件.

    由于条件 (a.)(b.) 可以确定唯一的三次埃尔米特插值多项式（笔者只能从结果反推这个结论应该是正确的，严格证明找相应文献）. 所以有两种做法，一种是直接把条件代入三次多项式通式方程求解待定系数，还有一种是借鉴拉格朗日插值基底的性质.

    .. math::
        :nowrap:

        \begin{align}
            l_{t}(x)=\dfrac{\prod\limits_{i=1}^{n\setminus t}(x-x_{i})}{\prod\limits_{i=1}^{n\setminus t}(x_t-x_{i})}
            \Longrightarrow
            l_t(x_s)=
            \begin{cases}
            1,\ \ t=s
            \\
            0,\ \ t\neq s
            \end{cases}
            \notag
        \end{align}

    基底满足条件：

        .. math:: \alpha_0(x_0)=1,\ \alpha_0(x_1)=0,\ \alpha_0'(x_0)=\alpha_0'(x_1)=0

        .. math:: \alpha_1(x_0)=0,\ \alpha_1(x_1)=1,\ \alpha_1'(x_0)=\alpha_1'(x_1)=0

        .. math:: \beta_0(x_0)=\beta_0(x_1)=0,\ \beta_0'(x_0)=1,\ \beta_0'(x_1)=0

        .. math:: \beta_1(x_0)=\beta_1(x_1)=0,\ \beta_1'(x_0)=0,\ \beta_1'(x_1)=1

    接下来借鉴性（很有难度）地构造出基底范式：
        .. math:: \aleph(x)=(ax+b)\times\left(\dfrac{x-x_i}{x_t-x_i}\right)^2

    通过范式模板，用待定系数把基底满足条件代入后可解出基底：
        .. math:: \alpha_0(x)=\left(1+2\times\dfrac{x-x_0}{x_1-x_0}\right)\left(\dfrac{x-x_1}{x_0-x_1}\right)^2

        .. math:: \alpha_1(x)=\left(1+2\times\dfrac{x-x_1}{x_0-x_1}\right)\left(\dfrac{x-x_0}{x_1-x_0}\right)^2

        .. math:: \beta_0(x)=(x-x_0)\left(\dfrac{x-x_1}{x_0-x_1}\right)^2 \ \ \ \ \beta_1(x)=(x-x_1)\left(\dfrac{x-x_0}{x_1-x_0}\right)^2

    余项：
        .. math:: R(x)=\dfrac{f^{(4)}(\xi)}{4!}(x-x_0)^2(x-x_1)^2


|

================
Ⅳ. 代码
================

.. code-block:: python
    :caption: 三点三次埃尔米特插值.py
    :linenos:

    '''
    # System --> Linux & Python3.8.0
    # File ----> 三点三次埃尔米特插值.py
    # Author --> Illusionna
    # Create --> 2024/2/20 21:18:17
    '''
    # -*- Encoding: UTF-8 -*-


    class HERMITE:
        def __init__(self, X:list, Y:list, middlePointDerivative:float) -> None:
            self.__X = X
            self.__Y = Y
            self.__middlePointDerivative = middlePointDerivative
            self.__pos = False

        def Coefficients(self) -> list:
            """公有函数: 计算埃尔米特多项式系数."""
            L = [self.__Y[0]]
            tmpA = (self.__Y[1] - self.__Y[0]) / (self.__X[1] - self.__X[0])
            L.append(tmpA)
            tmpB = ((self.__Y[2] - self.__Y[1]) / (self.__X[2] - self.__X[1]) - tmpA) / (self.__X[2] - self.__X[0])
            L.append(tmpB)
            C = (self.__middlePointDerivative - tmpA - tmpB*(self.__X[1] - self.__X[0])) / ((self.__X[1] - self.__X[0]) * (self.__X[1] - self.__X[2]))
            L.append(C)
            self.__coefficients = L
            self.__pos = True
            return L

        def Interpolate(self, x:float) -> float:
            """公有函数: 插值."""
            if self.__pos == True:
                result = 1*self.__coefficients[0]
                tmp = 1
                for index in range(0, len(self.__X), 1):
                    value = x - self.__X[index]
                    tmp = tmp * value
                    result = result + tmp * self.__coefficients[index+1]
                return result
            else:
                self.Coefficients()
                return self.Interpolate(x)


    if __name__ == '__main__':
        print('\033[H\033[J')

        X = [1/4, 1, 9/4]
        Y = [1/8, 1, 27/8]
        middlePointDerivative = 1.5     # 中间点导数, 即 x1 处的导数.

        obj = HERMITE(X, Y, middlePointDerivative)
        print(f'埃尔米特多项式系数:\n{obj.Coefficients()}')
        print(f'x = 1.6 处的插值结果: {obj.Interpolate(1.6)}')



.. code-block:: python
    :caption: 两点三次埃尔米特插值.py
    :linenos:

    '''
    # System --> Linux & Python3.8.0
    # File ----> 两点三次埃尔米特插值.py
    # Author --> Illusionna
    # Create --> 2024/2/20 21:46:23
    '''
    # -*- Encoding: UTF-8 -*-


    class HERMITE:
        def __init__(
                self,
                X:list,
                Y:list,
                derivative:list
        ) -> None:
            self.__X = X
            self.__Y = Y
            self.__derivative = derivative

        def Interpolate(self, x:float) -> 'function':
            """
            公有函数: 返回插值函数地址.
            """
            (x0, x1) = (self.__X[0], self.__X[1])
            (y0, y1) = (self.__Y[0], self.__Y[1])
            (diff0, diff1) = (self.__derivative[0], self.__derivative[1])
            alpha0 = lambda x: (1 + 2*((x-x0) / (x1-x0))) * ((x-x1) / (x0-x1))**2
            alpha1 = lambda x: (1 + 2*((x-x1) / (x0-x1))) * ((x-x0) / (x1-x0))**2
            beta0 = lambda x: (x-x0) * ((x-x1) / (x0-x1))**2
            beta1 = lambda x: (x-x1) * ((x-x0) / (x1-x0))**2
            H = y0*alpha0(x) + y1*alpha1(x) + diff0*beta0(x) + diff1*beta1(x)
            return H


    if __name__ == '__main__':
        print('\033[H\033[J')

        obj = HERMITE(
            X = [0, 1],
            Y = [0, 1],
            derivative = [-1, -4]
        )

        f:'function' = lambda x: obj.Interpolate(x)

        nodes:list = [-0.25, 0.25, 0.75, 1.25]
        predictions:list = list(map(f, nodes))

        print(f'节点: {nodes}\n\n插值: {predictions}')