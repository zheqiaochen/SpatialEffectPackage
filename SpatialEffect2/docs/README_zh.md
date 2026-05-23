# SpatialEffect2（中文说明）

`SpatialEffect2` 用于估计空间干预下、按距离索引的
**平均边际化效应（AME）**，适用于存在未知干扰（unknown
interference）的设置。

主要功能： - 支持点干预与多边形干预 - 支持 `edge` / `disk` / `donut`
三种距离定义 - 支持随机化分位数区间（`Per.CI`） - 支持 Conley 空间 HAC
区间（`Conley.CI`） - 支持局部多项式平滑与偏差修正

教程已直接写到 pkgdown 首页： - [pkgdown/index.md](pkgdown/index.md)

在包根目录生成站点：

``` r
pkgdown::build_site()
```

英文主文档： - [README.md](README.md)

------------------------------------------------------------------------

## 安装

本地源码安装：

``` r
install.packages("/path/to/SpatialEffect2", repos = NULL, type = "source")
```

或在包根目录执行：

``` bash
R CMD INSTALL SpatialEffect2
```

------------------------------------------------------------------------

## 快速接口

核心估计函数：

``` r
SpatialEffect(...)
```

距离区间检验：

``` r
SpatialEffectTest(result.list, test.range, smooth = 0, alpha = 0.05)
```

S3 方法：

``` r
summary(result, dVec.range = c(1, 5))
plot(result, smooth = TRUE, ci.type = "both")
```

------------------------------------------------------------------------

## 输入格式模板

### 1）点干预（`ras_Z = NULL`）

``` r
Zdata <- data.frame(
  id = 1:6,
  x = c(10, 15, 23, 30, 41, 44),
  y = c(12, 18, 20, 28, 33, 40),
  treat = c(1, 0, 1, 0, 1, 0)
)
```

调用时需要指定：

``` r
treatment = "treat",
x_coord_Z = "x",
y_coord_Z = "y"
```

### 2）多边形干预（`ras_Z` 为 `sf`）

``` r
result <- SpatialEffect(
  ras = ras,
  ras_Z = intervention_sf,
  Zdata = data.frame(treat = c(1, 0, 1, 0)),
  treatment = "treat",
  dVec = seq(0, 20, by = 1),
  cType = "donut"
)
```

`SpatialEffect2` 会自动用多边形质心作为节点坐标。

### 3）无栅格时的 kriging 模式（`ras = NULL`）

``` r
Ydata <- data.frame(
  x = runif(200, 0, 100),
  y = runif(200, 0, 100),
  outcome = rnorm(200)
)
```

配套参数：

``` r
ras = NULL,
Ydata = Ydata,
outcome = "outcome",
x_coord_Y = "x",
y_coord_Y = "y"
```

------------------------------------------------------------------------

## 参数说明（`SpatialEffect`）

### 必填参数

- `Zdata`：每行一个干预节点。
- `treatment`：`Zdata` 中二元处理变量列名。
- `dVec`：估计 AME 的距离网格。

此外： - 点干预必须给 `x_coord_Z`, `y_coord_Z`。 - 多边形干预给 `ras_Z`
即可（坐标自动推导）。

### 数据与几何

- `ras`：[`terra::SpatRaster`](https://rdrr.io/pkg/terra/man/SpatRaster-class.html)
  结果面；若为 `NULL` 则走 kriging。
- `Ydata`, `outcome`, `x_coord_Y`, `y_coord_Y`：kriging 输入。
- `ras_Z`：多边形干预（`sf` / `SpatRaster` / `SpatVector`）。

### 距离定义

- `cType`：
  - `"edge"`：只看边界
  - `"disk"`：看累计 `[0,d]`
  - `"donut"`：相邻累计环差分（常用）
- `dist.metric`：`"Euclidean"` 或 `"Geodesic"`。
- `numpts`, `evalpts`, `only.unique`：圈层采样控制。

### 设计校正与加权

- `covs`：协变量矩阵（`nz x p` 或 `nz*length(dVec) x p`）。
- `prob_treatment`：倾向得分列名（IPW）。
- `blockvar`, `clustvar`：置换推断时分层/聚类变量。

### 推断参数

- `per.se`：是否计算随机化区间 `Per.CI`。
- `conley.se`：是否计算 Conley 区间 `Conley.CI`。
- `kernel`：`"uni"`, `"tri"`, `"epa"`。
- `cutoff`：Conley 核带宽。
- `alpha`：显著性水平。
- `edf`：是否启用有效自由度调整。
- `nPerms`：置换次数。

### 平滑参数

- `smooth`：是否平滑。
- `bw`, `bw_debias`：平滑与偏差修正带宽（`NULL` 时自动 CV）。
- `bias_correction`：是否偏差修正。
- `smooth.conley.se`：是否计算平滑曲线标准误。
- `conf.band`：是否计算平滑曲线统一置信带。

### 计算参数

- `m`, `lambda`：kriging 参数。
- `n_threads`：并行线程数。

------------------------------------------------------------------------

## 输出格式

`SpatialEffect(...)` 返回 `"SpatialEffect"` 对象，常见元素：

- `AME_est`：`data.frame`，列为 `d`, `taud_est`
- `Per.CI`：`length(dVec) x 2`
- `Conley.SE`：长度为 `length(dVec)` 的向量
- `Conley.CI`：`length(dVec) x 2`
- `AME_est_smoothed`：平滑结果（`d`, `tau_smoothed`）
- `smoothed.Conley.CI` / `smoothed.Conley.CB`
- `Parameters`：内部参数与中间量

`SpatialEffectTest(...)` 输出： - `test.stat`：给定区间 AME 求和统计量 -
`test.CI`：sharp null 下置换区间

------------------------------------------------------------------------

## 实务流程建议

1.  先做几何体检。

- 保证 CRS 一致。
- 保证干预和结果面的空间重叠。

2.  先跑保守基线。

- 建议从 `cType = "donut"`、`smooth = 0`、`per.se = 1`、`conley.se = 1`
  开始。

3.  再加平滑。

- 先看未平滑结果稳定性，再启用 `smooth = 1`。
- 若自动带宽不稳，固定 `bw` / `bw_debias`。

4.  面向政策区间做检验。

- 用 `SpatialEffectTest(result, test.range = c(a, b))`
  直接评估关键距离段。

5.  做稳健性检查。

- 对比 `edge` / `disk` / `donut`。
- 改 `cutoff`、`kernel`、`dVec` 步长，检查结论是否稳定。

------------------------------------------------------------------------

## 常见问题

### 1）结果大量 `NA`

- 检查 CRS。
- 检查空间重叠。
- 检查 `Zdata` 与干预几何行数和顺序。
- 对 `donut` 适当调大 `dVec` 步长。

### 2）平滑不稳定或小样本矩阵奇异

- 增大 `bw` / `bw_debias`。
- 降低协变量共线性。
- 放宽距离网格密度。

### 3）运行慢

- 增大 `n_threads`。
- 降低 `nPerms`。
- 降低圈采样强度（`numpts`, `evalpts`）。

------------------------------------------------------------------------

## 参考文献

Wang, Y., Samii, C., Chang, H., & Aronow, P. M. (2024).  
Design-Based Inference for Spatial Experiments under Unknown
Interference.
