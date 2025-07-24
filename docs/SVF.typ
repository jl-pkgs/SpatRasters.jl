#import "@local/modern-cug-report:0.1.2": *
// #import "@preview/modern-cug-report:0.1.1": *
#counter(heading).update(2)
#let delta(x) = $Delta #x$

#show: (doc) => template(doc, 
  footer: "CUG水文气象学2024",
  header: "蒸散发的基本原理")


- O：人所在点为，坐标原点

- M：山顶所在点，$(  )$

- $S$: 人所在点的坡面角度

- $A$：坡面方位角

// $arrow("OM")$与坡面的方位角为$alpha$，

== 1 坡面法向量

正东为x轴，正北为y轴，z轴垂直于地面
$
arrow(n) = (sin S cos A, sin S sin A,cos S)  
$
- $S$: 坡度， 坡度法向量与z轴的夹角为天顶角，天顶角与坡度相等
- $A$: 数学上的方位角，逆时针为正，正东为0

== 2 光线方向
$
arrow(r)=( sin theta cos Phi, sin theta sin Phi, cos theta) 
$

- $theta$ : 太阳光线天顶角
- $Phi$   : 太阳光线方位角 （同A，数学上的方位角）

== 3 光线与人所在坡面法向量夹角 $alpha$
$
  cos(alpha) = arrow(n) arrow(r) 
            = sin S sin theta cos( Phi - A)+cos S cos theta 
$

立体角面积在某一平面的投影 =  $d omega * cos alpha$（图#[@fig_1]c） ，$d omega = sin(theta) d Phi d theta$。


#figure(
  image("../Figures/Figure1_示意图.png", width: 40%),
  caption: []
) <fig_1>

#pagebreak()

== 4 

$
  cos(alpha) d omega &= [sin(S) sin(theta) cos( Phi - A)+cos(S) cos(theta)] sin(theta) d Phi d theta \ 
  &= [sin(S) sin(theta)^2 cos( Phi - A) + cos(S) sin(theta) cos(theta)] d Phi d theta
$

- 常量：坡面方位角$A$、坡面坡度$S$

- 变量：太阳光线的天顶角$theta in [0, H]$，$Phi in [0, 2pi]$。其中H为地形最高点对应的天顶角，$H = pi/2 - "坡度"$。注意$H$与人所在点的坡面天顶角不同，需要根据山顶的高度进行指定。

对$ integral_0^(H) cos(alpha) d omega$进行积分在$theta in [0, H]$范围内进行积分：

#mitex(`$$
  \int_0^{H} \left[\sin(S)\sin^2(\theta)\cos(\Phi - A) + \cos(S)\sin(\theta)\cos(\theta)\right] d\theta \\
= 
\frac{1}{2} \sin(S)\cos(\Phi - A) \left( H - \frac{\sin(2H)}{2} \right) + \frac{1}{2} \cos(S)\sin^2 H
$$`)

#v(-0.7em)
#beamer-block[详细推导过程见：https://chatgpt.com/c/6881b3b4-f82c-8012-b2a6-ef910d1e9675。上式只有H是变量。]

因此，天空可视因子：

$
  "SVF" = 1/pi integral_0^(2pi) integral_0^(H) cos(alpha) d omega \
    = 1/ (2pi) integral_0^(2pi) sin(S) cos(Phi - A)[H - sin(H) cos(H)] + cos(S)sin(H)^2 d Phi
    // = integral_0^(2pi) integral_0^(pi/2 - s) [sin(S) sin(theta)^2 cos( Phi - A) + cos(S) sin(theta) cos(theta)]  d theta d Phi
$
