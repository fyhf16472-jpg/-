# 圆锥测距非线性最小二乘偏差验证

本仓库用于验证圆锥测距模型（如 GNSS-A 场景）下非线性最小二乘解算的二阶偏差，并提供可重复的 MATLAB/Octave 仿真代码。

## 仿真代码
- `matlab/cone_range_bias_sim.m`：无外部依赖的 Gauss-Newton 求解与蒙特卡洛偏差验证脚本，附带二阶 Hessian 近似的理论偏差计算。

运行示例（MATLAB 或 GNU Octave）：

```matlab
results = cone_range_bias_sim(struct('sigma',0.5,'nTrials',5000,'symmetric',false));
```

## 文档
- `docs/数值验证圆锥测距偏差.md`：对原稿的精简修订版，包含模型说明、偏差近似公式以及复现实验的操作步骤。
