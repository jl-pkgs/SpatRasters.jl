// #import "@preview/modern-cug-report:0.1.1": *
#import "@local/modern-cug-report:0.1.2": *

// #counter(heading).update(0)
#let out = "out"
#let delta(x) = $Delta #x$

#show: doc => template(
  doc,
  footer: "CUG水文气象学2025",
  header: "空间插值",
)

= 1 *不规则网格的双线性插值*

#beamer-block[#par(leading: 0.9em)[
参考`pyresample`的实现方法：\ https://github.com/pytroll/pyresample/blob/main/pyresample/bilinear/_base.py]]

#h(2em)
对于不规则四边形的四个角点 $P_1, P_2, P_3, P_4$，目标点 $(out_x, out_y)$ 可以表示为：

$
  P(s,t) = (1-s)(1-t) P_1 + s(1-t)P_2 + (1-s) t P_3 + s t P_4
$

其中 s 和 t 是需要求解的参数，范围都
在 [0,1] 区间内。

目标点 ($"out"_x$, $"out"_y$) 满足：

$
  "out"_x & = (1-s)(1-t) x₁ + s(1-t) x₂ + (1-s)t x₃ + s t x₄ \
          & = x₁ - s x₁ - t x₁ + s t x₁ + s x₂ - s t x₂ + t x₃ - s t x₃ + s t x₄ \
          & = x₁ + s(x₂ - x₁) + t(x₃ - x₁) + s t (x₁ - x₂ - x₃ + x₄)
$

$
  "out"_y & = (1-s)(1-t) y₁ + s(1-t) y₂ + (1-s)t y₃ + s t y₄ \
          & = y₁ + s(y₂ - y₁) + t(y₃ - y₁) + s t (y₁ - y₂ - y₃ + y₄)
$

定义：
$
  x₂₁ = x₂ - x₁, y₂₁ = y₂ - y₁ \
  x₃₁ = x₃ - x₁, y₃₁ = y₃ - y₁ \
  x₄₂ = x₄ - x₂, y₄₂ = y₄ - y₂
$

$
  "out"_x - x₁ = s x₂₁ + t x₃₁ + s t(x₄₂ - x₃₁) \
  // "out"_x - x₁ - t x₃₁ = s[x₂₁ + t(x₄₂ - x₃₁)] \
  s = ("out"_x - x₁ - t x₃₁) / [x₂₁ + t(x₄₂ - x₃₁)]
$ <eq_s>

将$s$代入$"out"_y$：

$
  "out"_y - y₁ = s y₂₁ + t y₃₁ + s t(y₄₂ - y₃₁) \
  "out"_y - y₁ = [("out"_x - x₁ - t x₃₁) / (x₂₁ + t(x₄₂ - x₃₁))] y₂₁ + t y₃₁
  + t[("out"_x - x₁ - t x₃₁) / (x₂₁ + t(x₄₂ - x₃₁))] (y₄₂ - y₃₁) \
$

两边同乘以分母：
$
  ("out"_y - y₁)[x₂₁ + t(x₄₂ - x₃₁)] = \
  ("out"_x - x₁ - t x₃₁)y₂₁ + t y₃₁[x₂₁ + t(x₄₂ - x₃₁)]
  + t("out"_x - x₁ - t x₃₁)(y₄₂ - y₃₁)
$


#box-blue[
  $
    (out_y - y₁)x₂₁ + t(out_y - y₁)(x₄₂ - x₃₁) = \
    (out_x - x₁)y₂₁ - t x₃₁y₂₁ + t y₃₁x₂₁ + t²y₃₁(x₄₂ - x₃₁)
    + t(out_x - x₁)(y₄₂ - y₃₁) - t²x₃₁(y₄₂ - y₃₁)
  $
]

将上式整理成$a t^2 + b t + c = 0$的形式：

$
  [x₃₁(y₄₂ - y₃₁) - y₃₁(x₄₂ - x₃₁)] t^2 + \
  [(out_y - y₁) (x₄₂ - x₃₁) + x₃₁ y₂₁ - y₃₁ x₂₁ - (out_x - x₁) (y₄₂ - y₃₁)] t + \
  (out_y - y₁)x₂₁ - (out_x - x₁)y₂₁ = 0
$

#box-red[
  $
    a = x₃₁(y₄₂ - y₃₁) - y₃₁(x₄₂ - x₃₁) = x₃₁y₄₂ - y₃₁x₄₂ \
    b = (out_y - y₁) (x₄₂ - x₃₁) + x₃₁ y₂₁ - y₃₁ x₂₁ - (out_x - x₁) (y₄₂ - y₃₁) \
    quad = out_y (x₄₂ - x₃₁) - out_x (y₄₂ - y₃₁) - x₄₂ y₁ + x₃₁ y₂ + x₁ y₄₂ - y₃₁ x₂ \
    c = (out_y - y₁)x₂₁ - (out_x - x₁)y₂₁ = out_y x₂₁ - out_x y₂₁ + x₁y₂ - y₁x₂
  $
]

之后，将t代入式#[@eq_s]可得`s`。
