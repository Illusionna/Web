高次埃尔米特插值示例
======================

.. role:: blue
    :class: blue

.. role:: red
    :class: red

.. raw:: html

    <style>

        .blue {
            color:blue;
        }
        .red {
            color:red;
        }

    </style>


|


Ⅰ. 前言 --> Ⅱ. 示例


|

================
Ⅰ. 前言
================

涉及到牛顿插值均差系数，建议先阅读：

:doc:`../牛顿插值/Newton`

|

================
Ⅱ. 示例
================

计算已知点列和条件构成的插值多项式：

点列：

.. math:: \bm{\vec{x}}=(-1,\ 0,\ 1)\ \ \ \ \bm{\vec{y}}=(0,\ -4,\ -2)

条件：

.. math::
    :nowrap:

    \begin{align}
        \begin{cases}
        \bm{\vec{y'}}=({\rm Unknown},\ 0,\ 5)
        \\
        \bm{\vec{y''}}=({\rm Unknown},\ 6,\ {\rm Unknown})
        \end{cases}
        \notag
    \end{align}


在 :math:`x_1=0` 处出现了一阶导数 :math:`0` 和二阶导数 :math:`6`，共计 2 次，在 :math:`x_2=1` 处出现一阶导数 :math:`5`，共计 1 次.

所以构造如下均差表：

.. list-table:: 
  :widths: 2 2 6 6 5 5 5
  :header-rows: 1

  * - x
    - y
    - 一阶差商
    - 二阶差商
    - 三阶差商
    - 四阶差商
    - 五阶差商
  * - -1
    - :red:`0`
    -
    -
    -
    -
    -
  * - 0
    - :math:`-4`
    - (-4-0)/[0-(-1)]= :red:`-4`
    -
    -
    -
    -
  * - :blue:`0`
    - :blue:`-4`
    - [-4-(-4)]/[(0-0) :math:`\times1!`]= :math:`y_1'\div1!` =0
    - [0-(-4)]/[0-(-1)]= :red:`4`
    -
    -
    -
  * - :blue:`0`
    - :blue:`-4`
    - [-4-(-4)]/[(0-0) :math:`\times1!`]= :math:`y_1'\div1!` =0
    - (0-0)/[(0-0) :math:`\times2!`]= :math:`y_1''\div2!` =3
    - (3-4)/[0-(-1)]= :red:`-1`
    - 
    -
  * - 1
    - :math:`-2`
    - [-2-(-4)]/(1-0)=2
    - (2-0)/(1-0)=2
    - (2-3)/(1-0)=-1
    - [-1-(-1)]/[1-(-1)]= :red:`0`
    -
  * - :blue:`1`
    - :blue:`-2`
    - [-2-(-2)]/[(1-1) :math:`\times1!`]= :math:`y_2'\div1!` =5
    - (5-2)/(1-0)=3
    - (3-2)/(1-0)=1
    - [1-(-1)]/(1-0)=2
    - (2-0)/[1-(-1)]= :red:`1`

所以五次埃尔米特插值多项式为：

.. math:: H(x)=0\times1-4[x-(-1)]+4[x-(-1)](x-0)-[x-(-1)](x-0)^2+0+[x-(-1)](x-0)^3(x-1)

|